#!/usr/bin/env python

import sys
import argparse
import numpy as np
import csv
import pysam
from time import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from contextlib import redirect_stdout
from typing import Collection, Any

from dark.fasta import FastaReads
from dark.sam import samfile, samReferenceLengths
from dark.utils import baseCountsToStr, pct

BASES = "ACGT"
TRANSITIONS = "AG GA CT TC".split()
TRANSVERSIONS = "AT TA AC CA GT TG GC CG".split()


def makeParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Print SAM/BAM per-site nucleotide composition and coverage statistics. "
            "Note that no output is produced for sites with no read coverage."
        ),
    )

    parser.add_argument(
        "samfile",
        help="The SAM/BAM file to read.",
    )

    parser.add_argument(
        "--referenceId",
        metavar="ID",
        help=(
            "A reference sequence id. If not given and the SAM file has "
            "just one reference, that one will be used. Otherwise, a "
            "reference id must be given."
        ),
    )

    parser.add_argument(
        "--fasta",
        help="Optionally give a FASTA file in which the reference can be found.",
    )

    parser.add_argument(
        "--window",
        "-n",
        type=int,
        default=250,
        help=(
            "The number of genome sites to have each subprocess handle. Large values "
            "will have a major performance impact when --maxDepth is also high (e.g., "
            ">10,000). This seems to be due to an unknown pysam/htslib slowdown (based "
            "on samtools v1.21, as of 2025-01-31)."
        ),
    )

    parser.add_argument(
        "--workers",
        type=int,
        help=(
            "The number of subprocesses to run in parallel. If not specified, the "
            "process pool will typically use all available cores, but see "
            "https://docs.python.org/3/library/concurrent.futures.html"
            "#processpoolexecutor for details."
        ),
    )

    parser.add_argument(
        "--noOffsets",
        action="store_false",
        dest="printOffsets",
        help="Do not print per-site details of base counts.",
    )

    parser.add_argument(
        "--noStats",
        action="store_false",
        dest="printStats",
        help="Do not print final average, standard deviation etc. statistics.",
    )

    parser.add_argument(
        "--statsFile",
        help=(
            "The file to print overall statistics to. If not specified, standard "
            "output will be used."
        ),
    )

    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help=(
            "Print progress information (to standard error). 0 = don't print progress "
            "information. 1 = print information after each genome region is processed. "
            "2 = also print information during the processing of each region."
        ),
    )

    parser.add_argument(
        "--minDepth",
        type=int,
        help=(
            "The minimum read depth to report on. Sites with a lower depth will "
            "not be reported."
        ),
    )

    parser.add_argument(
        "--maxDepth",
        type=int,
        default=1e6,
        help="The maximum read depth to request from pysam for each site.",
    )

    parser.add_argument(
        "--minSiteMutationRate",
        type=float,
        help=(
            "The minimum site mutation rate to report on. Sites with a lower mutation "
            "rate will not be reported. The mutation rate is the fraction of bases at "
            "the site that are not the commonest base."
        ),
    )

    parser.add_argument(
        "--noFilter",
        action="store_false",
        dest="filter",
        help="This behaviour is now the default and this option is ignored.",
    )

    parser.add_argument(
        "--tsv",
        action="store_true",
        help=(
            "Produce TSV output (use --sep , if you want CSV). "
            "Note that overall SAM file statistics will be printed at the bottom of "
            "the output unless you also specify --noStats or --statsFile."
        ),
    )

    parser.add_argument(
        "--sep",
        default="\t",
        help="The separator to use when --tsv is used. Ignored otherwise.",
    )

    return parser


def getReferenceId(
    filename: str, references: Collection[str], referenceId: str | None
) -> str:
    """
    Get the reference id.  This is either the reference id from the command line
    (if given and if present in the SAM file) or the single reference from the
    SAM file (if the SAM file contains just one reference).
    """
    if referenceId:
        if referenceId not in references:
            raise ValueError(
                f"Reference id {referenceId} not found in {filename!r} "
                f"Known references are {', '.join(sorted(references))}."
            )
    else:
        if len(references) == 1:
            referenceId = list(references)[0]
        else:
            raise ValueError(
                f"Multiple references found in {filename!r}: "
                f"{', '.join(sorted(references))}. Use --referenceId "
                "to specify the one you want."
            )

    return referenceId


def getReferenceSeq(fasta: str | None, referenceId: str) -> str | None:
    """
    Look for the reference id in the FASTA file and return its sequence if found.
    If no FASTA is given, return None.
    """
    if fasta:
        for read in FastaReads(fasta):
            if read.id == referenceId:
                return read.sequence
        else:
            raise ValueError(
                f"Reference id {referenceId!r} was not found in reference FASTA "
                f"file {fasta!r}."
            )
    else:
        return None


def processSubsection(
    filename: str,
    start: int,
    window: int,
    minDepth: int | None,
    maxDepth: int,
    minSiteMutationRate: float | None,
    collectOffsets: bool,
    verbosity: int,
    referenceId: str,
    referenceSeq: str | None,
    returnList: bool,
    startTime: float,
) -> dict[str, Any]:
    """
    Process a subsection of the genome.

    @param filename: The SAM/BAM file to read.
    @param start: The start offset of the genome window we process.
    @param window: The size of the genome section we should process.
    @param miDepth: The minimum read depth to report on. Sites with lower depth will
        be ignored.
    @param maxDepth: The maximum read depth to consider.
    @param minSiteMutationRate: The minimal mutation rate (i.e., fraction) required at
        sites in order that they are reported. Sites with lower mutation fractions will
        be ignored.
    @param collectOffsets: Whether to collect genome offset information.
    @param verbosity: The verbosity level.
    @param referenceId: The name of the reference sequence.
    @param referenceSeq: The genome of the reference (if known).
    @param returnList: If True, the rows we return will be printed (by our caller) as
        CSV, so we return lists. Else, we return a string for each row of output.
    @param startTime: The time overall processing was launched.
    """
    readIds = set()
    rows = []
    totalBases = totalReads = totalMutationCount = 0
    totalMutationPairCounts = defaultdict(int)
    depths = []
    startTime = time()

    with samfile(filename) as sam:
        for column in sam.pileup(
            start=start,
            end=start + window,
            reference=referenceId,
            truncate=True,
            max_depth=maxDepth,
        ):
            assert isinstance(column, pysam.PileupColumn)
            columnStartTime = time()
            refOffset = column.reference_pos
            assert start <= refOffset < start + window

            bases = dict.fromkeys(BASES, 0)
            nReads = 0
            for read in column.pileups:
                if read.query_position is not None:
                    nReads += 1
                    readIds.add(read.alignment.query_name)
                    assert read.alignment.query_sequence is not None
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base] += 1

            if minDepth is not None and nReads < minDepth:
                if verbosity > 1:
                    print(
                        f"Site {refOffset + 1}: skipped due to too few reads ({nReads:,}).",
                        file=sys.stderr,
                    )
                continue

            totalBases += nReads
            mutationCount = 0
            refBase = referenceSeq[refOffset] if referenceSeq else None

            def sortKey(base):
                return (bases[base], int(base == refBase))

            # Find the most-common base (with a preference for the reference
            # base (if any) in case of draws) and count mutations as the
            # total number of bases that are not the commonest.
            commonestBase = max(bases, key=sortKey)

            mutationRate = 1.0 - bases[commonestBase] / nReads
            if minSiteMutationRate is not None and mutationRate < minSiteMutationRate:
                if verbosity > 1:
                    print(
                        f"Site {refOffset + 1}: skipped due to low mutation rate "
                        f"({mutationRate:.3f}).",
                        file=sys.stderr,
                    )
                continue

            mutationPairCounts = defaultdict(int)
            for base in BASES:
                if base != commonestBase and (count := bases[base]):
                    mutationPairCounts[commonestBase + base] += count
                    mutationCount += count

            depths.append(nReads)

            totalMutationCount += mutationCount
            for pair, count in mutationPairCounts.items():
                totalMutationPairCounts[pair] += count

            if collectOffsets:
                if returnList:
                    row: list[int | str | float] = [refOffset + 1]
                    if refBase:
                        row.append(refBase)
                    row.extend(
                        [nReads]
                        + [bases[base] for base in BASES]
                        + [mutationCount, round(mutationRate, 6)]
                    )

                    totalTransitions = 0
                    for pair in TRANSITIONS:
                        count = mutationPairCounts.get(pair, 0)
                        row.append(count)
                        totalTransitions += count

                    totalTransversions = 0
                    for pair in TRANSVERSIONS:
                        count = mutationPairCounts.get(pair, 0)
                        row.append(count)
                        totalTransversions += count

                    row.extend((totalTransitions, totalTransversions))

                    rows.append(row)
                else:
                    out = f"{refOffset + 1} {nReads} {baseCountsToStr(bases)}"
                    if refBase:
                        out += f" Ref:{refBase}"
                    rows.append(out)

            if verbosity > 1:
                totalReads += nReads
                stopTime = time()
                elapsed = stopTime - columnStartTime
                totalElapsed = stopTime - startTime

                print(
                    f"Site {refOffset + 1}: {nReads:,} reads processed "
                    f"at {nReads / elapsed / 1e6:.2f} M read/s. "
                    f"Total bases: {totalBases:,} "
                    f"Overall: {totalReads / totalElapsed / 1e6:.2f} M read/s",
                    file=sys.stderr,
                )

    return {
        "depths": depths,
        "mutationPairCounts": totalMutationPairCounts,
        "readIds": readIds,
        "rows": rows,
        "startSite": start,
        "startTime": startTime,
        "stopTime": time(),
        "totalBases": totalBases,
        "totalMutationCount": totalMutationCount,
    }


def main() -> None:
    args = makeParser().parse_args()

    if not (args.printOffsets or args.printStats):
        sys.exit("You used both --noOffsets and --noStats, so there is no output!")

    referenceLens = samReferenceLengths(args.samfile)
    referenceId = getReferenceId(args.samfile, set(referenceLens), args.referenceId)
    referenceLen = referenceLens[referenceId]
    referenceSeq = getReferenceSeq(args.fasta, referenceId)
    depths = []
    printTSV = csv.writer(sys.stdout, delimiter=args.sep).writerow if args.tsv else None
    printRow = printTSV or print

    if printTSV:
        header = ["Site"]
        if referenceSeq:
            header.append("Reference")
        header.append("Depth")
        header.extend(BASES)
        header.append("Mutations")
        header.append("Variability")
        header.extend(f"{pair[0]}->{pair[1]}" for pair in TRANSITIONS + TRANSVERSIONS)
        header.extend(("Transitions", "Transversions"))
        printTSV(header)

    readIds = set()
    totalBases = totalMutationCount = 0
    # I don't use a collections.Counter in the following because they are very slow.
    mutationPairCounts = defaultdict(int)
    start = time()
    totalRunTime = maxRunTime = 0.0
    startSites = list(range(0, referenceLen, args.window))

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for result in executor.map(
            processSubsection,
            repeat(args.samfile),
            startSites,
            repeat(args.window),
            repeat(args.minDepth),
            repeat(args.maxDepth),
            repeat(args.minSiteMutationRate),
            repeat(args.printOffsets),
            repeat(args.verbosity),
            repeat(referenceId),
            repeat(referenceSeq),
            repeat(bool(printTSV)),
            repeat(start),
        ):
            waited = result["startTime"] - start
            runTime = result["stopTime"] - result["startTime"]
            totalRunTime += runTime
            if runTime > maxRunTime:
                maxRunTime = runTime

            if args.verbosity > 0:
                from_ = result["startSite"]
                to = min(from_ + args.window, referenceLen)
                print(
                    f"Genome region {from_ + 1}-{to} processed in {runTime:.2f} "
                    f"seconds (waited {waited:.2f} seconds to start).",
                    file=sys.stderr,
                )

            if args.printOffsets:
                for row in result["rows"]:
                    printRow(row)

            readIds.update(result["readIds"])
            totalBases += result["totalBases"]
            depths.extend(result["depths"])
            for pair, count in result["mutationPairCounts"].items():
                mutationPairCounts[pair] += count
            totalMutationCount += result["totalMutationCount"]

    if args.verbosity > 0:
        stop = time()
        try:
            speedup = totalRunTime / (stop - start)
        except ZeroDivisionError:
            speedup = 0.0

        print(
            f"Finished in {int(stop - start)} seconds. Total serial processing would "
            f"have been {int(totalRunTime)} seconds. The slowest process took "
            f"{maxRunTime:.2f} seconds. Speed up = {speedup:.2f}",
            file=sys.stderr,
        )

    if args.printStats:
        out = open(args.statsFile, "w") if args.statsFile else sys.stdout
        with redirect_stdout(out):
            print("SAM file:", args.samfile)
            print("Max depth considered:", args.maxDepth)
            print("Min depth required:", args.minDepth)
            print("Min (site) mutation rate:", args.minSiteMutationRate)
            print("Reference id:", referenceId)
            print("Reference length:", referenceLen)
            print("Number of reads found:", len(readIds))
            print("Reference genome coverage:", pct(len(depths), referenceLen))
            print("Coverage depth:")
            print("  Min:", min(depths, default=0))
            print("  Max:", max(depths, default=0))

            if depths:
                print(f"  Mean: {np.mean(depths):.3f}")
                print(f"  Median: {np.median(depths):.3f}")
                print(f"  SD: {np.std(depths):.3f}")

            omr = totalMutationCount / totalBases if totalBases else 0.0
            print(f"Overall mutation rate: {omr:.6f}")

            totalTransitions = totalTransversions = 0
            print("Sorted mutation counts (* = transversion):")
            countWidth = None

            def key(pair):
                return mutationPairCounts[pair], pair

            for pair in sorted(mutationPairCounts, key=key, reverse=True):
                count = mutationPairCounts[pair]
                if countWidth is None:
                    countWidth = len(str(count))
                if pair in TRANSITIONS:
                    totalTransitions += count
                    what = ""
                else:
                    assert pair in TRANSVERSIONS
                    totalTransversions += count
                    what = " *"
                print(f"  {pair[0]} -> {pair[1]}: {count:{countWidth}d}{what}")

            print("Total transitions:", totalTransitions)
            print("Total transversions:", totalTransversions)

        out.close()


if __name__ == "__main__":
    main()

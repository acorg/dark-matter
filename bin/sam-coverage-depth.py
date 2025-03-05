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
            "A reference sequence id. If not given and the SAM/BAM file has just one "
            "reference, that one will be used. Otherwise, a reference id must be given."
        ),
    )

    parser.add_argument(
        "--fasta",
        help=(
            "Optionally give a FASTA file in which the reference can be found. "
            "If you do not provide a reference sequence, mutations will be counted "
            "as differences from the most-common nucleotide at each genome site."
        ),
    )

    parser.add_argument(
        "--diffsFrom",
        choices=("reference", "commonest"),
        default="reference",
        help=(
            "The base from which to assess differences. If 'reference', differences "
            "from the reference base (at each genome site) will be treated as "
            "mutations. If 'commonest', differences from the most common base will be "
            "treated as mutations (with ties broken in favour of the reference base, "
            "if a reference is given). It is important that you understand the "
            "difference between these two and choose the right option."
        ),
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
            "the site that are not the reference or commonest base (depending on "
            "the --diffsFrom option)."
        ),
    )

    parser.add_argument(
        "--maxSiteMutationRate",
        type=float,
        help=(
            "The maximum site mutation rate to report on. Sites with a higher mutation "
            "rate will not be reported. The mutation rate is the fraction of bases at "
            "the site that are not the reference or commonest base (depending on "
            "the --diffsFrom option)."
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
            "Produce TSV output (use --sep , if you want CSV). Note that overall "
            "SAM/BAM file statistics will be printed at the bottom of the output "
            "unless you also specify --noStats or --statsFile. Columns will hopefully "
            "be self-explanatory, but a few clarifying comments will probably be "
            "helpful. 1) The 'Base' column indicates the nucleotide base that was used "
            "to count mutations. The base used will depend on the value of --diffsFrom. "
            "See the help output for that option. 2) The 'Reference' column will have "
            "an empty value if no reference sequence was provided. 3) The 'Site' column "
            "gives 1-based indexes into the reference genome."
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
    (if given and if present in the SAM/BAM file) or the single reference from the
    SAM/BAM file (if it contains just one reference).

    @param filename: The SAM/BAM filename.
    @param references: The names of all references in the SAM/BAM file.
    @param referenceId: The id of the wanted reference, if provided on the command line.
    @return: The reference id.
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


def getReferenceSeq(filename: str, referenceId: str) -> str:
    """
    Look for the reference id in the FASTA file and return its sequence if found.
    If no FASTA is given, return None.

    @param filename: The FASTA filename.
    @param referenceId: The ID of the wanted reference.
    @return: The reference nucleotide sequence.
    """
    for read in FastaReads(filename):
        if read.id == referenceId:
            return read.sequence
    else:
        raise ValueError(
            f"Reference id {referenceId!r} was not found in reference FASTA "
            f"file {filename!r}."
        )


def processSubsection(
    filename: str,
    start: int,
    window: int,
    minDepth: int | None,
    maxDepth: int,
    minSiteMutationRate: float | None,
    maxSiteMutationRate: float | None,
    collectOffsets: bool,
    verbosity: int,
    referenceId: str,
    referenceSeq: str | None,
    returnList: bool,
    startTime: float,
    diffsFrom: str,
) -> dict[str, Any]:
    """
    Process a subsection of the genome.

    @param filename: The SAM/BAM file to read.
    @param start: The start offset of the genome window we process.
    @param window: The size of the genome section we should process.
    @param miDepth: The minimum read depth to report on. Sites with lower depth will
        be ignored.
    @param maxDepth: The maximum read depth to consider.
    @param minSiteMutationRate: The minimum mutation rate (i.e., fraction) required at
        sites in order that they are reported. Sites with lower mutation fractions will
        be ignored.
    @param maxSiteMutationRate: The maximum mutation rate (i.e., fraction) allowed at
        sites in order that they are reported. Sites with higher mutation fractions will
        be ignored.
    @param collectOffsets: Whether to collect genome offset information.
    @param verbosity: The verbosity level.
    @param referenceId: The name of the reference sequence.
    @param referenceSeq: The genome of the reference (if known).
    @param returnList: If True, the rows we return will be printed (by our caller) as
        CSV, so we return lists. Else, we return a string for each row of output.
    @param startTime: The time overall processing was launched.
    @param diffsFrom: The base from which to assess differences. If 'reference',
        differences from the reference base (at each genome site) will be treated as
        mutations. If 'commonest', differences from the most common base will be
        treated as mutations (with ties broken in favour of the reference base,
        if a reference is given).
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

            if diffsFrom == "reference":
                assert refBase is not None
                fromBase = refBase
            else:
                def sortKey(base: str) -> tuple[int, int]:
                    """
                    Return a key that can be used to find the base with the highest count
                    (via max). If there is no reference, the second value in the returned
                    tuple will always be False (because refBase will be None). As a result,
                    there is no mechanism for deterministically resolving ties (i.e.,
                    equally-common nucleotide counts at a genome site. If there is a
                    reference, ties are broken in favour of the reference base.
                    """
                    return (bases[base], int(base == refBase))

                fromBase = max(bases, key=sortKey)

            mutationRate = 1.0 - bases[fromBase] / nReads

            if minSiteMutationRate is not None and mutationRate < minSiteMutationRate:
                if verbosity > 1:
                    print(
                        f"Site {refOffset + 1}: skipped due to low mutation rate "
                        f"({mutationRate:.5f}).",
                        file=sys.stderr,
                    )
                continue

            if maxSiteMutationRate is not None and mutationRate > maxSiteMutationRate:
                if verbosity > 1:
                    print(
                        f"Site {refOffset + 1}: skipped due to high mutation rate "
                        f"({mutationRate:.5f}).",
                        file=sys.stderr,
                    )
                continue

            mutationPairCounts = defaultdict(int)
            for base in BASES:
                if base != fromBase and (count := bases[base]):
                    mutationPairCounts[fromBase + base] += count
                    mutationCount += count

            depths.append(nReads)

            totalMutationCount += mutationCount
            for pair, count in mutationPairCounts.items():
                totalMutationPairCounts[pair] += count

            if collectOffsets:
                if returnList:
                    row: list[int | str | float | None] = (
                        [refOffset + 1, refBase, nReads]
                        # Don't use bases.values() in the following, because we
                        # need bases in ACGT order.
                        + [bases[base] for base in BASES]
                        + [fromBase, mutationCount, round(mutationRate, 6)]
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
                    if referenceSeq:
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
    referenceSeq = getReferenceSeq(args.fasta, referenceId) if args.fasta else None

    if args.diffsFrom == "reference" and referenceSeq is None:
        sys.exit(
            "If you use --diffsFrom 'reference' (the default), you must use --fasta to "
            "give a FASTA file containing the reference sequence. Else use --diffsFrom "
            "'commonest' to calculate mutations based on the most-common nucleotide "
            "at each genome site."
        )

    depths = []
    printTSV = csv.writer(sys.stdout, delimiter=args.sep).writerow if args.tsv else None
    printRow = printTSV or print

    if printTSV:
        printTSV(
            ("Site", "Reference", "Depth") +
            tuple(BASES) +
            ("Base", "Mutations", "Variability") +
            tuple(f"{pair[0]}->{pair[1]}" for pair in TRANSITIONS + TRANSVERSIONS) +
            ("Transitions", "Transversions"),
        )

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
            repeat(args.maxSiteMutationRate),
            repeat(args.printOffsets),
            repeat(args.verbosity),
            repeat(referenceId),
            repeat(referenceSeq),
            repeat(bool(printTSV)),
            repeat(start),
            repeat(args.diffsFrom),
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
            print("SAM/BAM file:", args.samfile)
            print("Mutations are relative to:", args.diffsFrom)
            print("Max depth considered:", args.maxDepth)
            print("Min depth required:", args.minDepth)
            print("Min (site) mutation rate:", args.minSiteMutationRate)
            print("Max (site) mutation rate:", args.maxSiteMutationRate)
            print("Reference id:", referenceId)
            print("Reference length:", referenceLen)
            print("Reference sequence given:", "yes" if referenceSeq else "no")
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

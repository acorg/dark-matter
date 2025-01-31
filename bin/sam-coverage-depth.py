#!/usr/bin/env python

import sys
import argparse
import numpy as np
import csv
from time import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from contextlib import redirect_stdout

from dark.fasta import FastaReads
from dark.sam import samfile, samReferenceLengths
from dark.utils import baseCountsToStr, pct

BASES = "ACGT"


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
        "--maxDepth",
        type=int,
        default=1e6,
        help="The maximum read depth to request from pysam for each site.",
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


def getReferenceId(filename: str, references: set[str], referenceId: str | None) -> str:
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
    args: argparse.Namespace,
    start: int,
    referenceId: str,
    referenceSeq: str | None,
    returnCSV: bool,
    startTime: float,
):
    """
    Process a subsection of the genome.
    """
    readIds = set()
    rows = []
    totalBases = totalReads = totalMutationCount = 0
    mutationPairCounts = defaultdict(int)
    baseCounts = []
    startTime = time()

    with samfile(args.samfile) as sam:
        for column in sam.pileup(
            start=start,
            end=start + args.window,
            reference=referenceId,
            truncate=True,
            max_depth=args.maxDepth,
        ):
            columnStartTime = time()
            refOffset = column.reference_pos
            assert start <= refOffset < start + args.window
            pileups = column.pileups
            nReads = column.nsegments

            bases = dict.fromkeys(BASES, 0)
            for read in pileups:
                if read.query_position is not None:
                    readIds.add(read.alignment.query_name)
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base] += 1

            baseCount = sum(bases.values())
            baseCounts.append(baseCount)
            totalBases += baseCount

            if referenceSeq:
                refBase = referenceSeq[refOffset]

                def sortKey(base):
                    return (bases[base], int(base == refBase))

                # Find the most-common base (with a preference for the
                # reference base in case of draws) and count mutations as the
                # total number of bases that are not the commonest.
                commonestBase = sorted(bases, key=sortKey, reverse=True)[0]
                mutationCount = 0
                for base in BASES:
                    if base != commonestBase and (baseCount := bases[base]):
                        mutationPairCounts[commonestBase + base] += baseCount
                        mutationCount += baseCount
                totalMutationCount += mutationCount
            else:
                refBase = None

            if args.printOffsets:
                if returnCSV:
                    row = [refOffset + 1]
                    if refBase:
                        row.extend((refBase, mutationCount))
                    row.extend(bases[base] for base in BASES)
                    rows.append(row)
                else:
                    out = f"{refOffset + 1} {baseCount} " f"{baseCountsToStr(bases)}"
                    if refBase:
                        out += f" Ref:{refBase}"
                    rows.append(out)

            if args.verbosity > 1:
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
        "baseCounts": baseCounts,
        "mutationPairCounts": mutationPairCounts,
        "readIds": readIds,
        "rows": rows,
        "startSite": start,
        "startTime": startTime,
        "stopTime": time(),
        "totalBases": totalBases,
        "totalMutationCount": totalMutationCount,
    }


def main():
    args = makeParser().parse_args()

    if not (args.printOffsets or args.printStats):
        exit("You used both --noOffsets and --noStats, so there is no output!")

    referenceLens = samReferenceLengths(args.samfile)
    referenceId = getReferenceId(args.samfile, set(referenceLens), args.referenceId)
    referenceLen = referenceLens[referenceId]
    referenceSeq = getReferenceSeq(args.fasta, referenceId)
    baseCounts = []
    writerow = csv.writer(sys.stdout, delimiter=args.sep).writerow if args.tsv else None

    if writerow:
        row = ["Site"]
        if referenceSeq:
            row.extend(("Reference", "Mutations"))
        writerow(row + list(BASES))

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
            repeat(args),
            startSites,
            repeat(referenceId),
            repeat(referenceSeq),
            repeat(bool(writerow)),
            repeat(start),
        ):
            waited = result["startTime"] - start
            runTime = result["stopTime"] - result["startTime"]
            totalRunTime += runTime
            if runTime > maxRunTime:
                maxRunTime = runTime

            if args.verbosity > 0:
                from_ = result['startSite']
                to = min(from_ + args.window, referenceLen)
                print(
                    f"Genome region {from_ + 1}-{to} processed in {runTime:.2f} "
                    f"seconds, after waiting {waited:.2f} seconds to start.",
                    file=sys.stderr,
                )

            func = writerow or print
            for row in result["rows"]:
                func(row)

            readIds.update(result["readIds"])
            totalBases += result["totalBases"]
            baseCounts.extend(result["baseCounts"])
            for pair, count in result["mutationPairCounts"].items():
                mutationPairCounts[pair] += count
            totalMutationCount += result["totalMutationCount"]

    stop = time()
    try:
        speedup = totalRunTime / (stop - start)
    except ZeroDivisionError:
        speedup = 0.0

    if args.verbosity > 0:
        print(
            f"Finised in {int(stop - start)} seconds. Total serial processing would "
            f"have been {int(totalRunTime)} seconds. The slowest process took "
            f"{maxRunTime:.2f} seconds. Speed up = {speedup:.2f}",
            file=sys.stderr,
        )

    if args.printStats:
        out = open(args.statsFile, "w") if args.statsFile else sys.stdout
        with redirect_stdout(out):
            referenceLength = referenceLens[referenceId]
            print("SAM file: %s" % args.samfile)
            print("Max depth: %d" % args.maxDepth)
            print("Reference id: %s" % referenceId)
            print("Reference length: %d" % referenceLength)
            print("Number of reads found: %d" % len(readIds))
            print("Bases covered: %s" % pct(len(baseCounts), referenceLength))
            print("Min coverage depth: %d" % min(baseCounts, default=0))
            print("Max coverage depth: %d" % max(baseCounts, default=0))

            if baseCounts:
                print("Mean coverage depth: %.3f" % np.mean(baseCounts))
                print("Median coverage depth: %.3f" % np.median(baseCounts))
                print("Coverage depth s.d.: %.3f" % np.std(baseCounts))
                if referenceSeq:
                    print(f"Mutation rate: {totalMutationCount / totalBases:.6f}")

                    def key(pair):
                        return mutationPairCounts[pair], pair

                    print("Mutation counts:")
                    for pair in sorted(mutationPairCounts, key=key, reverse=True):
                        ref, mut = pair
                        print(f"  {ref} -> {mut}: {mutationPairCounts[pair]}")

        out.close()


if __name__ == "__main__":
    main()

#!/usr/bin/env python

import sys
import argparse
import numpy as np
import csv
from time import time
from collections import defaultdict
import concurrent.futures
from itertools import repeat

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
        help="The number of genome sites to have each subprocess handle.",
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
        help="Do not print final average and standard deviation statistics.",
    )

    parser.add_argument(
        "--statsFile",
        help=(
            "The file to print overall statistics to. If not specified, standard "
            "output will be used."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress information (to standard error).",
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
            print(
                f"Reference id {referenceId!r} was not found in reference FASTA file {fasta!r}",
                file=sys.stderr,
            )
            sys.exit(1)
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

            if args.verbose:
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
        "startSite": start,
        "startTime": startTime,
        "stopTime": time(),
        "rows": rows,
        "totalBases": totalBases,
        # "totalReads": totalReads,
        "baseCounts": baseCounts,
        "totalMutationCount": totalMutationCount,
        "mutationPairCounts": mutationPairCounts,
        "readIds": readIds,
    }


def main():
    parser = makeParser()
    args = parser.parse_args()

    if not (args.printOffsets or args.printStats):
        print(
            "You have used both --noOffsets and --noStats, so there is no output!",
            file=sys.stderr,
        )
        sys.exit(1)

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
    mutationPairCounts = defaultdict(int)
    start = time()
    totalRunTime = maxRunTime = 0.0
    startSites = list(range(0, referenceLen, args.window))

    with concurrent.futures.ProcessPoolExecutor() as executor:
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

            print(
                f"Window {result['startSite']} produced in {runTime:.2f} seconds, after waiting "
                f"{waited:.2f} seconds to start.",
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

    print(
        f"Finised in {int(stop - start)} seconds. Total serial processing would "
        f"have been {int(totalRunTime)} seconds. The slowest process took "
        f"{maxRunTime:.2f} seconds. Speed up = {speedup:.2f}",
        file=sys.stderr,
    )

    if args.printStats:
        out = open(args.statsFile, "w") if args.statsFile else sys.stdout
        referenceLength = referenceLens[referenceId]
        print("SAM file: %s" % args.samfile, file=out)
        print("Max depth: %d" % args.maxDepth, file=out)
        print("Reference id: %s" % referenceId, file=out)
        print("Reference length: %d" % referenceLength, file=out)
        print("Number of reads found: %d" % len(readIds), file=out)
        print("Bases covered: %s" % pct(len(baseCounts), referenceLength), file=out)
        print("Min coverage depth: %d" % min(baseCounts, default=0), file=out)
        print("Max coverage depth: %d" % max(baseCounts, default=0), file=out)

        if baseCounts:
            print("Mean coverage depth: %.3f" % np.mean(baseCounts), file=out)
            print("Median coverage depth: %.3f" % np.median(baseCounts), file=out)
            print("Coverage depth s.d.: %.3f" % np.std(baseCounts), file=out)
            if referenceSeq:
                print(f"Mutation rate: {totalMutationCount / totalBases:.6f}", file=out)

                def key(pair):
                    return mutationPairCounts[pair], pair

                print("Mutation counts:", file=out)
                for pair in sorted(mutationPairCounts, key=key, reverse=True):
                    ref, mut = pair
                    print(f"  {ref} -> {mut}: {mutationPairCounts[pair]}", file=out)

        out.close()


if __name__ == "__main__":
    main()

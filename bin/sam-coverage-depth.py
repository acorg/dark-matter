#!/usr/bin/env python

import sys
import argparse
import numpy as np
import csv
from time import time
from collections import defaultdict

from dark.fasta import FastaReads
from dark.sam import samfile, SAMFilter, samReferences, UnknownReference
from dark.utils import baseCountsToStr, pct


def makeParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Print SAM/BAM per-site nucleotide composition and coverage statistics. "
            "Note that no output is produced for sites with no read coverage."
        ),
    )

    SAMFilter.addFilteringOptions(parser, samfileIsPositional=True)

    parser.add_argument(
        "--fasta",
        help="Optionally give a FASTA file in which the reference can be found.",
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
        )
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
        help="This behaviour is now the default and this option is ignored."
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


def main():
    BASES = "ACGT"
    parser = makeParser()
    args = parser.parse_args()

    if not (args.printOffsets or args.printStats):
        print(
            "You have used both --noOffsets and --noStats, so there is no output!",
            file=sys.stderr,
        )
        sys.exit(1)

    samFilter = SAMFilter.parseFilteringOptions(args)

    if samFilter.referenceIds and len(samFilter.referenceIds) > 1:
        print(
            "Only one reference id can be given. To calculate coverage for more "
            "than one reference, run this script multiple times.",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        referenceLengths = samFilter.referenceLengths()
    except UnknownReference:
        referenceId = samFilter.referenceIds.pop()
        referenceIds = samReferences(args.samfile)
        print(
            "Reference %r does not appear in SAM file %s. Known "
            "references are: %s."
            % (referenceId, args.samfile, ", ".join(sorted(referenceIds))),
            file=sys.stderr,
        )
        sys.exit(1)

    baseCounts = []
    writerow = csv.writer(sys.stdout, delimiter=args.sep).writerow if args.tsv else None

    with samfile(args.samfile) as sam:
        if samFilter.referenceIds:
            # No need to check if the given reference id is in referenceLengths
            # because the samFilter.referenceLengths call above catches that.
            referenceId = samFilter.referenceIds.pop()
        else:
            if len(referenceLengths) == 1:
                referenceId = list(referenceLengths)[0]
            else:
                print(
                    "SAM file %r contains %d references (%s). Only one "
                    "reference id can be analyzed at a time. Please use "
                    "--referenceId to specify the one you want examined."
                    % (
                        args.samfile,
                        len(referenceLengths),
                        ", ".join(sorted(referenceLengths)),
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)

        if args.fasta:
            for read in FastaReads(args.fasta):
                if read.id == referenceId:
                    referenceSeq = read.sequence
                    break
            else:
                print(
                    f"Reference id {referenceId!r} was not found in reference FASTA file {args.fasta!r}",
                    file=sys.stderr,
                )
                sys.exit(1)
        else:
            referenceSeq = None

        if writerow:
            row = ["Site"]
            if referenceSeq:
                row.extend(("Reference", "Mutations"))
            writerow(row + list(BASES))

        readIds = set()
        totalBases = totalReads = totalMutationCount = 0
        mutationPairCounts = defaultdict(int)
        startTime = time()

        for column in sam.pileup(
            reference=referenceId,
            truncate=True,
            max_depth=args.maxDepth,
        ):
            columnStartTime = time()
            refOffset = column.reference_pos
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
                if writerow:
                    row = [refOffset + 1]
                    if refBase:
                        row.extend((refBase, mutationCount))
                    row.extend(bases[base] for base in BASES)
                    writerow(row)
                else:
                    out = f"{refOffset + 1} {baseCount} " f"{baseCountsToStr(bases)}"
                    if refBase:
                        out += f" Ref:{refBase}"
                    print(out)

            if args.verbose:
                totalReads += nReads
                stop = time()
                elapsed = stop - columnStartTime
                totalElapsed = stop - startTime

                print(
                    f"Site {refOffset + 1}: {nReads:,} reads processed "
                    f"at {nReads / elapsed / 1e6:.2f} M read/s. "
                    f"Total bases: {totalBases:,} "
                    f"Overall: {totalReads / totalElapsed / 1e6:.2f} M read/s",
                    file=sys.stderr,
                )

    if args.printStats:
        out = open(args.statsFile, "w") if args.statsFile else sys.stdout
        referenceLength = referenceLengths[referenceId]
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

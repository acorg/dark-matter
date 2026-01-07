#!/usr/bin/env python

import argparse
import re
import sys
from collections import defaultdict

import numpy as np
import polars as pl

from dark.filter import (
    addFASTAFilteringCommandLineOptions,
    parseFASTAFilteringCommandLineOptions,
)
from dark.reads import Reads
from dark.sam import SAMFilter, samfile


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Print SAM/BAM file coverage statistics.",
    )

    parser.add_argument(
        "--noFilter",
        action="store_true",
        dest="filter",
        default=False,
        help=(
            "Do not use our SAM filtering. This now the default since filtering "
            "can be extremely slow. Use --filter to turn filtering on. Note that "
            "if you use this option, any filtering option (other than --referenceId) "
            "you also specify that is provided by the SAMFilter.addFilteringOptions "
            "will be silently ignored!"
        ),
    )

    parser.add_argument(
        "--includeZeroes",
        action="store_true",
        help="Also write output for references that have no covering reads.",
    )

    parser.add_argument(
        "--filter",
        action="store_true",
        help="Use our SAM filtering. Note that this can be extremely slow.",
    )

    parser.add_argument(
        "--regex",
        action="store_true",
        help="Interpret the --referenceId arguments as regular expressions.",
    )

    parser.add_argument(
        "--match",
        action="store_true",
        help=(
            "Use 'match' when checking if references match regular expressions. "
            "This will cause the regular expression to be checked against the start "
            "of the reference name. If not given, the regular expression must "
            "match the entire reference name not just its prefix."
        ),
    )

    parser.add_argument(
        "--listReferences",
        action="store_true",
        help=(
            "Print the list of matching reference IDs and then exit. "
            "This makes it easy to adjust reference regular expressions "
            "to see what work would be done without doing it."
        ),
    )

    parser.add_argument(
        "--minReferenceLength",
        type=int,
        help=(
            "The minimum length reference to use. This can be useful if you give a "
            "reference regular expression that matches too many references and you "
            "want to also limit them by length. Or you can just give a minimal length "
            "and that way pick up all chromosomes without needing to name them."
        ),
    )

    parser.add_argument(
        "--tsv",
        help="The name of a TSV file to save results to.",
    )

    parser.add_argument(
        "--excel",
        help="The name of an Excel file to save results to.",
    )

    addFASTAFilteringCommandLineOptions(parser)
    SAMFilter.addFilteringOptions(parser, samfileIsPositional=True)

    return parser.parse_args()


def main() -> None:
    args = get_args()

    # We don't have a file of reads, we just want a read filter that we can use
    # to filter the SAM file query sequences and to get reference lengths from.
    reads = parseFASTAFilteringCommandLineOptions(args, Reads())
    samFilter = SAMFilter.parseFilteringOptions(args, reads.filterRead)

    if args.filter:
        filtering = True

        def filterRead(read):
            return not (read.is_del or read.is_refskip) and samFilter.filterAlignment(
                read.alignment
            )
    else:
        filtering = False

        def filterRead(read):
            return True

    referenceLengths = samFilter.referenceLengths(getAllReferences=True)

    if args.referenceId:
        if args.regex:
            regex = re.compile("|".join(map(re.escape, args.referenceId)))
            matcher = regex.match if args.match else regex.fullmatch
            wantedIds = set(id_ for id_ in referenceLengths if matcher(id_))

            if not wantedIds:
                sys.exit(
                    "None of the reference ids in the BAM file match the reference "
                    "id regular expression(s) you specified!"
                )

        else:
            wantedIds = set(args.referenceId)

            if not (wantedIds & set(referenceLengths)):
                sys.exit("None of the given reference ids occur in the input BAM file!")

    else:
        wantedIds = set(referenceLengths)

    if args.minReferenceLength is not None:
        n_before = len(wantedIds)
        wantedIds = set(
            id_ for id_ in wantedIds if referenceLengths[id_] >= args.minReferenceLength
        )

        if not wantedIds:
            sys.exit(
                f"None of the {n_before} potentially-included reference(s) has length "
                f"at least {args.minReferenceLength}."
            )

    if args.listReferences:
        print("\n".join(sorted(wantedIds)))
        sys.exit(0)

    coveredOffsetCount: dict[str, dict[int, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    coveringReads = defaultdict(set)

    with samfile(args.samfile) as sam:
        for column in sam.pileup():
            referenceOffset = column.reference_pos
            for read in column.pileups:
                referenceId = read.alignment.reference_name
                if referenceId in wantedIds and (not filtering or filterRead(read)):
                    coveredOffsetCount[referenceId][referenceOffset] += 1
                    coveringReads[referenceId].add(read.alignment.query_name)

    rows = []

    for referenceId in sorted(wantedIds):
        offsetsCovered = len(coveredOffsetCount[referenceId])

        if offsetsCovered == 0 and not args.includeZeroes:
            continue

        referenceLength = referenceLengths[referenceId]

        nonZero = list(coveredOffsetCount[referenceId].values())
        assert len(nonZero) <= referenceLength
        coverageCounts = np.zeros(referenceLength, dtype=int)
        coverageCounts[: len(nonZero)] = nonZero

        rows.append(
            (
                referenceId,
                referenceLength,
                len(coveringReads[referenceId]),
                offsetsCovered,
                offsetsCovered / referenceLength,
                np.mean(coverageCounts),
                np.median(coverageCounts),
                np.min(coverageCounts),
                np.max(coverageCounts),
            )
        )

    if args.tsv or args.excel:
        schema = (
            "reference id",
            "reference length",
            "covering reads",
            "offsets covered",
            "fraction covered",
            "mean depth",
            "median depth",
            "min depth",
            "max depth",
        )

        df = pl.DataFrame(rows, schema=schema, orient="row")

        if args.tsv:
            df.write_csv(args.tsv, separator="\t")

        if args.excel:
            df.write_excel(args.excel)
    else:
        for row in rows:
            (
                referenceId,
                referenceLength,
                coveringReads,
                offsetsCovered,
                fraction,
                meanDepth,
                medianDepth,
                minDepth,
                maxDepth,
            ) = row

            print(
                f"{referenceId}: length {referenceLength}, covering "
                f"reads {coveringReads}, covered sites {offsetsCovered} "
                f"({fraction * 100.0:.4f}%), mean coverage depth {meanDepth:.4f}, "
                f"median coverage depth {medianDepth:.4f} "
                f"(min: {minDepth}, max: {maxDepth})"
            )


if __name__ == "__main__":
    main()

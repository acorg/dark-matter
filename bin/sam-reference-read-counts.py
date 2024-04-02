#!/usr/bin/env python

import sys
import argparse

from dark.sam import ReadsSummary


def main(args):
    """
    Print SAM/BAM file reference read counts.

    @param args: An argparse namespace with information about parsed
        command-line options.
    """
    if args.topReferenceIdsFile and args.sortBy not in ("count", "coverage"):
        print(
            "--topReferenceIdsFile only makes sense when sorting by count or coverage",
            file=sys.stderr,
        )
        sys.exit(1)

    summary = ReadsSummary(args.samFile)

    print(
        summary.toString(
            excludeZeroes=args.excludeZeroes,
            excludeIfNoAdditional=args.excludeIfNoAdditional,
            sortBy=args.sortBy,
        )
    )

    # Write out the sorted read ids of the top reference.
    if args.topReferenceIdsFile:
        with open(args.topReferenceIdsFile, "w") as fp:
            print("\n".join(sorted(summary.sortedReferences[0].readIds)), file=fp)


def makeParser():
    """
    Make an argument parser.

    @return: An C{argparse} argument parser.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Print SAM/BAM file reference sequence and unmapped read counts.",
    )

    parser.add_argument(
        "samFile", metavar="FILENAME", help="The name of a SAM/BAM alignment file."
    )

    parser.add_argument(
        "--sortBy",
        choices=("count", "coverage", "name"),
        default="count",
        help=(
            "Use 'count' to sort the output by decreasing read count (i.e., "
            "the reference with the highest number of matching reads is "
            "printed first), or 'coverage' to sort by decreasing genome "
            "coverage, or 'name' to sort by reference name."
        ),
    )

    parser.add_argument(
        "--excludeZeroes",
        action="store_true",
        help="Do not print references unmatched by any read.",
    )

    parser.add_argument(
        "--excludeIfNoAdditional",
        action="store_true",
        help=(
            "Do not print references that are mapped by no additional reads "
            "(i.e., reads other than those matching already-printed "
            "references). This is useful when --sortBy count is used, "
            "making it possible to not print references only matched by "
            "reads matching references that have already been printed."
        ),
    )

    parser.add_argument(
        "--topReferenceIdsFile",
        metavar="FILE",
        help=(
            "The file to write the (sorted) ids of the reads for the "
            "reference with the highest number of matching reads. Only "
            "valid when --sortBy count is used."
        ),
    )

    return parser


if __name__ == "__main__":
    main(makeParser().parse_args())

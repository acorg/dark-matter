#!/usr/bin/env python

import argparse
import csv
import sys

from summarize_ranges.ranges import RangeSummarizer

from dark.dna import AMBIGUOUS
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


def formatRanges(rs: RangeSummarizer, highAdjust: int) -> str:
    result = []
    for low, high in rs.ranges():
        high -= highAdjust
        if low == high:
            result.append(str(low))
        else:
            result.append(f"{low}-{high}")

    return ",".join(result)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin, write information about the "
            "counts of characters in the sequences to stdout."
        )
    )

    parser.add_argument(
        "--ignoreCase",
        "-i",
        action="store_true",
        help="Map all characters to uppercase.",
    )

    parser.add_argument(
        "--sep", default="\t", help="The output field separator character."
    )

    parser.add_argument(
        "--omitIds",
        action="store_true",
        help="Do not print sequence ids as the first field.",
    )

    parser.add_argument(
        "--omitHeader",
        dest="includeHeader",
        action="store_false",
        help="Do not print a header.",
    )

    parser.add_argument(
        "--ambiguous",
        action="store_true",
        help="If --chars is not given, print counts for all ambiguous DNA codes.",
    )

    parser.add_argument(
        "--omitZeroes",
        action="store_true",
        help="Do not print chars whose counts are zero.",
    )

    parser.add_argument(
        "--locations",
        action="store_true",
        help="Also print the 1-based genome locations of the matched characters.",
    )

    parser.add_argument(
        "--sort",
        action="store_true",
        help=(
            "Sort the characters for each sequence by decreasing count (with draws "
            "broken in favour of alphabetically increasing character)."
        ),
    )

    parser.add_argument(
        "--python",
        action="store_true",
        help=(
            "Print 0-based Python closed-open location ranges (i.e., the upper bound "
            "is not included) when --locations is used."
        ),
    )

    parser.add_argument(
        "--chars",
        "-c",
        help="The characters to count. Default is all IUPAC nucleotide codes.",
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()

    omitZeroes = args.omitZeroes
    addLocations = args.locations
    highAdjust = 0 if args.python else 1
    omitIds = args.omitIds
    reads = parseFASTACommandLineOptions(args)
    writerow = csv.writer(sys.stdout, delimiter=args.sep).writerow

    if args.chars:
        chars = args.chars
    else:
        if args.ambiguous:
            chars = "".join(
                sorted(c for (c, codons) in AMBIGUOUS.items() if len(codons) > 1)
            )
        else:
            chars = "".join(sorted(AMBIGUOUS))

    if args.ignoreCase:
        chars = chars.upper()

    if args.includeHeader:
        header = [] if omitIds else ["ID"]
        header.extend(("Char", "Count", "Length", "Fraction"))
        if addLocations:
            header.append("Location(s)")
        writerow(header)

    def key(row):
        return (row[1], row[0]) if omitIds else (row[2], row[1])

    for read in reads:
        length = len(read)
        sequence = read.sequence.upper() if args.ignoreCase else read.sequence
        rows = []
        for char in chars:
            if addLocations:
                locations = []
                rs = RangeSummarizer()

            count = 0
            for i, c in enumerate(sequence):
                if c == char:
                    count += 1
                    if addLocations:
                        rs.add(str(i + 1))
                        locations.append(i + 1)

            if count or not omitZeroes:
                row = [] if omitIds else [read.id]
                row.extend([char, count, length, f"{count / length:.4f}"])
                if addLocations:
                    row.append(formatRanges(rs, highAdjust))
                rows.append(row)

        if args.sort:
            rows.sort(key=key, reverse=True)

        for row in rows:
            writerow(row)


if __name__ == "__main__":
    main()

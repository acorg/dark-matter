#!/usr/bin/env python

import sys
import re
import csv
from typing import Union

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions

OVERLAPPING_INDICATOR = "+"

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin and a sequence pattern, write information "
            "about which sequences and where the pattern matches to stdout."
        )
    )

    parser.add_argument(
        "--regex",
        required=True,
        help=(
            "The regular expression pattern to match. For matching, input "
            "sequences will be converted to uppercase unless --noUpper is used."
        ),
    )

    parser.add_argument(
        "--noUpper",
        dest="upper",
        action="store_false",
        help="Do not convert sequences to uppercase for the matching.",
    )

    parser.add_argument(
        "--overlapping", action="store_true", help="Include overlapping matches."
    )

    parser.add_argument(
        "--noId",
        dest="printId",
        action="store_false",
        help="Do not print sequence ids.",
    )

    parser.add_argument(
        "--noEnd",
        dest="printEnd",
        action="store_false",
        help="Do not print the end offset of matches.",
    )

    parser.add_argument(
        "--noSequence",
        dest="printSequence",
        action="store_false",
        help="Do not print sequence match regions.",
    )

    parser.add_argument(
        "--noOverlapping",
        dest="printOverlapping",
        action="store_false",
        help=(
            f"Do not print a {OVERLAPPING_INDICATOR!r} indicator to "
            f"identify matches that overlap the previously printed one."
        ),
    )

    parser.add_argument(
        "--zeroBased",
        action="store_true",
        help="Print half-open zero-based offsets, as in Python.",
    )

    parser.add_argument(
        "--count",
        "-c",
        action="store_true",
        help="Just print the count of the number of matches, if the count is non-zero.",
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    regex = re.compile(args.regex)
    upper = args.upper
    printEnd = args.printEnd
    printId = args.printId
    printSequence = args.printSequence
    printOverlapping = args.printOverlapping
    overlapping = args.overlapping
    startInc = int(not args.zeroBased)
    writerow = csv.writer(sys.stdout, delimiter="\t").writerow

    for read in reads:
        sequence = read.sequence.upper() if upper else read.sequence
        startPos = 0
        matches = []
        while True:
            match = regex.search(sequence, pos=startPos)
            if match:
                start, end = match.start(), match.end()
                startPos = start + 1 if overlapping else end
                matches.append(match)
            else:
                break

        if args.count:
            if matches:
                if printId:
                    writerow((read.id, len(matches)))
                else:
                    print(len(matches))
        else:
            lastEnd = -1
            for match in matches:
                start, end = match.start(), match.end()
                row: list[Union[int, str]] = [start + startInc]
                if printEnd:
                    row.append(end)
                if printId:
                    row.append(read.id)
                if printSequence:
                    row.append(match.group())
                if printOverlapping and start < lastEnd:
                    row.append(OVERLAPPING_INDICATOR)
                writerow(row)
                lastEnd = end

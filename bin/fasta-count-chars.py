#!/usr/bin/env python

import sys
import csv

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin, write information about the "
            "number of various characters to stdout."
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
        help="Do not print sequence ids as the first field",
    )

    parser.add_argument("--chars", "-c", required=True, help="The characters to count.")

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    chars = args.chars.upper() if args.ignoreCase else args.chars
    reads = parseFASTACommandLineOptions(args)

    writerow = csv.writer(sys.stdout, delimiter=args.sep).writerow

    for read in reads:
        length = len(read)
        sequence = read.sequence.upper() if args.ignoreCase else read.sequence
        for char in chars:
            count = sequence.count(char)
            out = [] if args.omitIds else [read.id]
            writerow(out + [char, count, length, f"{count / length:.4f}"])

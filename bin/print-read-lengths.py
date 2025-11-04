#!/usr/bin/env python

"""
Read a FASTA or FASTQ file (or read stdin) and print a line for each sequence,
with the length of the sequence followed by its name.
"""

import argparse

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin, write the sequence lengths and ids to stdout."
        ),
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        print(len(read), read.id)


if __name__ == "__main__":
    main()

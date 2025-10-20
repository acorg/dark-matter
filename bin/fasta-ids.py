#!/usr/bin/env python

import argparse

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


def main():
    parser = argparse.ArgumentParser(
        description="Given FASTA on stdin, write the sequence ids to stdout.",
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        print(read.id)


if __name__ == "__main__":
    main()

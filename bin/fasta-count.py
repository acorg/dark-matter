#!/usr/bin/env python

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=("Given FASTA on stdin, print the number of sequences to stdout."),
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    count = 0
    for read in reads:
        count += 1

    print(count)

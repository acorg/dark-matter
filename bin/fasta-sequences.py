#!/usr/bin/env python

from __future__ import print_function

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Given FASTA on stdin, write the sequences to stdout.'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        print(read.sequence)

#!/usr/bin/env python

from __future__ import print_function

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Given FASTA on stdin, write the ids, sequences, and '
            'quality strings on a single TAB-separated line to stdout.'))

    parser.add_argument('--separator', default='\t',
                        help=('The character string to separate ids from '
                              'sequences and quality strings (if any)'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    sep = args.separator
    reads = parseFASTACommandLineOptions(args)

    # Duplicate code a little so as not to test for quality on each read.
    if args.fastq:
        for read in reads:
            print(read.id + sep + read.sequence + sep + read.quality)
    else:
        for read in reads:
            print(read.id + sep + read.sequence)

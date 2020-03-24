#!/usr/bin/env python

from __future__ import print_function

import hashlib

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Given FASTA on stdin, write the sequences to stdout, with '
            'various possibilities for formatting.'))

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '--md5', action='store_true', default=False,
        help='Print MD5 checksums of sequences.')

    group.add_argument(
        '--md5OneLine', action='store_true', default=False,
        help=('Print MD5 checksums of sequences, then the sequence id, then '
              'the sequence, then the quality string (if any) all '
              'TAB-separated.'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    # Below we duplicate the loops to avoid a test on every read.

    if args.md5:
        # Print just the MD5 sum of the sequences.
        md5 = hashlib.md5
        for read in reads:
            print(md5(read.sequence.encode('utf-8')).hexdigest())

    elif args.md5OneLine:
        # Print TAB-separated: MD5 sum of the sequences, read id, read
        # sequence, and sequence quality (if any).
        md5 = hashlib.md5
        if args.fasta:
            for read in reads:
                print('\t'.join(
                    md5(read.sequence.encode('utf-8')).hexdigest(),
                    read.id, read.sequence))
        else:
            for read in reads:
                print('\t'.join(
                    md5(read.sequence.encode('utf-8')).hexdigest(),
                    read.id, read.sequence, read.quality))
    else:
        for read in reads:
            print(read.sequence)

#!/usr/bin/env python

from __future__ import print_function

import hashlib

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin, write the sequences '
                     'to stdout.'))

    parser.add_argument(
        '--md5', action='store_true', default=False,
        help='Print MD5 checksums of seqeunces.')

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    # Duplicate the loop to avoid a test on every read.
    if args.md5:
        md5 = hashlib.md5
        for read in reads:
            print(md5(read.sequence.encode('utf-8')).hexdigest())
    else:
        for read in reads:
            print(read.sequence)

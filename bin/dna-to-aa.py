#!/usr/bin/env python

"""
Read DNA FASTA from stdin and print AA FASTA to stdout.  If a minimum
ORF length is given, only print AA sequences that have an ORF of at least
that length.
"""

import sys
import argparse

from dark.fasta import FastaReads


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert DNA to AA',
        epilog='Given DNA FASTA on stdin, output AA FASTA to stdout. '
        'Optionally, filter by minimum required ORF length.'
    )

    parser.add_argument(
        '--minORFLength', metavar='LEN', type=int, default=None,
        help='Translations to AA that do not contain an ORF of at least '
        'this length will not be produced.')

    args = parser.parse_args()
    reads = FastaReads(sys.stdin)
    write = sys.stdout.write
    minORFLength = args.minORFLength

    for read in reads:
        for aa in read.translations():
            if minORFLength is None or aa.maximumORFLength() >= minORFLength:
                write(aa.toString('fasta'))

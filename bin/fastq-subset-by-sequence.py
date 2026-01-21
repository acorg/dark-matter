#!/usr/bin/env python

"""
Given a pattern from sys.argv, read FASTQ from stdin,
and print FASTQ reads that contain an exact sequence match to stdout.
"""

from __future__ import print_function
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract a subset of FASTQ reads by sequence',
        epilog='Given a set of FASTQ sequence identifiers from sys.argv, '
        'read FASTQ from stdin, and print FASTQ that contain an '
        'exact sequence match to stdout.')

    parser.add_argument(
        '--pattern', nargs=1, 
         help='Wanted read sequence pattern. Only give one pattern.')

    args = parser.parse_args()
    pattern = Seq(args.pattern[0])

    found = []

    if pattern:
        for sequence in SeqIO.parse(sys.stdin, 'fastq'):
            if pattern in sequence.seq:
                found.append(sequence)

        if found:
            SeqIO.write(found, sys.stdout, 'fastq')
            print('Found %d sequences.' % len(found), file=sys.stderr)

        else:
            print('WARNING: Sequence with pattern %s not found' % (
                  pattern), file=sys.stderr)
    else:
        # No wanted patterns were given.
        parser.print_help(sys.stderr)
        sys.exit(1)
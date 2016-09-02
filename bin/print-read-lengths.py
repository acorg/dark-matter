#!/usr/bin/env python

"""
Read a FASTA or FASTQ file (or read stdin) and print a line for each sequence,
with the length of the sequence followed by its name.
"""

from __future__ import print_function

import sys
import argparse
from Bio import SeqIO


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Print lengths of reads in a FASTA or FASTQ file',
    )

    parser.add_argument(
        'reads', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
        help='The input sequence file (in either FASTA or FASTQ format). '
        'If not given, sequences will be read from stdin.')

    parser.add_argument(
        '--format', default='fasta', choices=['fasta', 'fastq'],
        help='The format of the input.')

    args = parser.parse_args()

    for seq in SeqIO.parse(args.reads, args.format):
        print('%d %s' % (len(seq), seq.description))

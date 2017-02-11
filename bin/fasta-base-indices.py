#!/usr/bin/env python

from __future__ import print_function

import os
import sys

from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.fastq import FastqReads


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Given FASTA on stdin, and a base to look for, write the 1-based '
            'indices of where the base occurs across *all* sequences to '
            'standard output. If standard output is a terminal these will be '
            'space separated, else newline separated.'),
        epilog=(
            'This can be used to find all columns of a FASTA multiple '
            'sequence alignment that contain gaps, or some other character.'))

    parser.add_argument(
        '--readClass', default='fasta', choices=('fasta', 'fastq', 'fasta-ss'),
        help='If specified, give the input FASTA type')

    parser.add_argument(
        '--base', default='-',
        help='The sequence base whose indices should be printed.')

    parser.add_argument(
        '--matchCase', default=False, action='store_true',
        help='If specified, sequence case will be considered.')

    args = parser.parse_args()

    if args.readClass == 'fastq':
        reads = FastqReads(sys.stdin)
    elif args.readClass == 'fasta':
        reads = FastaReads(sys.stdin, checkAlphabet=0)
    else:
        assert args.readClass == 'fasta-ss'
        reads = SSFastaReads(sys.stdin, checkAlphabet=0)

    result = set()
    target = args.base

    if args.matchCase:
        for read in reads:
            for index, base in enumerate(read.sequence, start=1):
                if base == target:
                    result.add(index)
    else:
        target = target.lower()
        for read in reads:
            for index, base in enumerate(read.sequence.lower(), start=1):
                if base == target:
                    result.add(index)

    separator = ' ' if os.isatty(1) else '\n'
    print(separator.join(map(str, sorted(result))))

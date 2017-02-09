#!/usr/bin/env python

from __future__ import print_function

import sys
from collections import defaultdict
from math import log10

from dark.aa import NAMES
from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.fastq import FastqReads
from dark.summarize import sequenceCategoryLengths


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin, write a summary of sequence base '
                     'categories to stdout. It is currently not possible to '
                     'specify the categories on the command line.'))

    parser.add_argument(
        '--baseType', default='nucl', choices=('nucl', 'prot'),
        help='The type of the bases in the input.')

    parser.add_argument(
        '--readClass', default='fasta', choices=('fasta', 'fastq', 'fasta-ss'),
        help='If specified, give the input FASTA type.')

    parser.add_argument(
        '--minLength', default=1, type=int,
        help=('If specified, stretches of reads that are less than this '
              'length will not be reported but will be summarized by an '
              'ellipsis.'))

    parser.add_argument(
        '--concise', action='store_true', default=False,
        help='If specified, do not show the individual sequence regions.')

    args = parser.parse_args()

    if args.readClass == 'fastq':
        # TODO: FastqReads should take a checkAlphabet argument, in the way
        # that FastaReads does.
        reads = FastqReads(sys.stdin)
    elif args.readClass == 'fasta':
        reads = FastaReads(sys.stdin, checkAlphabet=False)
    else:
        # args.readClass must be fasta-ss due to the 'choices' argument
        # passed to parser.add_argument value above.
        assert args.readClass == 'fasta-ss'
        reads = SSFastaReads(sys.stdin, checkAlphabet=False)

    if args.baseType == 'nucl':
        categories = {
            'A': 'nucl',
            'C': 'nucl',
            'G': 'nucl',
            'T': 'nucl',
            '-': 'gap',
        }
        default = 'ambiguous'
    else:
        categories = {}
        for name in NAMES:
            categories[name] = 'aa'
        categories['-'] = 'gap'
        default = 'ambiguous'

    categoryWidth = max(
        [len(category) for category in categories.values()] + [len(default)])

    minLength = args.minLength
    concise = args.concise

    for index, read in enumerate(reads, start=1):
        counts = defaultdict(int)
        readLen = len(read)
        width = int(log10(readLen)) + 1
        if not concise:
            summary = []
            append = summary.append
            offset = 1
        for (category, count) in sequenceCategoryLengths(
                read, categories, defaultCategory=default,
                minLength=minLength):
            counts[category] += count
            if not concise:
                append('    %*d %-*s (offset %*d)' %
                       (width, count, categoryWidth, category, width, offset))
                offset += count
        print('%d: %s (length %d)' % (index, read.id, readLen))
        for category in sorted(counts):
            count = counts[category]
            print('  %-*s: %*d (%6.2f%%)' %
                  (categoryWidth, category, width, count,
                   count / readLen * 100.0))
        if not concise:
            print('\n'.join(summary))

#!/usr/bin/env python

import sys
from Bio import SeqIO


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin and a minimum length ',
                     'write (FASTA) sequences of at least the given length '
                     'to stdout.'))

    parser.add_argument('len', type=int, help='The minimum sequence length')
    parser.add_argument('--verbose', type=bool, default=False,
                        help='If True, print information on found primers')

    args = parser.parse_args()
    reads = []
    short = count = 0

    for seqRecord in SeqIO.parse(sys.stdin, 'fasta'):
        count += 1
        if len(seqRecord) >= args.len:
            reads.append(seqRecord)
        else:
            short += 1

    if args.verbose:
        print >>sys.stderr, (
            'Read %d sequences. Sufficiently long: %d. Too short: %d.') % (
            count, count - short, short)

    SeqIO.write(reads, sys.stdout, 'fasta')

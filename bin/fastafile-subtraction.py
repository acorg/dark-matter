#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from dark import fastamanipulations

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate a fastafile, containing a specified '
        'set of reads',
        epilog='Given a list of readIds and a fastaFile, generate '
        'a new fastaFile based on these reads'
    )

    parser.add_argument(
        '--fasta', metavar='FASTA-files', type=list,
        help='a list of fastafiles which should be subtracted.')

    parser.add_argument(
        '--out', type=str, default=None,
        help='Where the output should be written to. If None, '
        'write to stdout.')

    args = parser.parse_args()

    reads = fastamanipulations.fastaSubtract(args.fasta)

    if not args.out:
        SeqIO.write(reads, sys.stdout, 'fasta')
    else:
        with open(args.out, 'w') as fp:
            SeqIO.write(reads, fp, 'fasta')

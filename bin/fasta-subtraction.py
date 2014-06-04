#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from dark import fastamanipulations



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate a fastafile, containing a specified set of reads',
        epilog='Given a list of readIds and a fastaFile, generate a new fastaFile '
        'based on these reads'
    )

    parser.add_argument(
        '--fasta1', metavar='FASTA-files', type=str,
        help='the FASTA file of sequences from which we wish to subtract.')

    parser.add_argument(
        '--readIds', metavar='readId-list',
        help='either a list with readIds or a file containing one readId per line.')

    parser.add_argument(
        '--present', default=True, type=bool,
        help='If True, all reads in readIds will be returned, else, all reads not in '
        'readIds will be returned.')

    parser.add_argument(
        '--out', type=str, default=None,
        help='Where the output should be written to. If None, write to stdout.')

    args = parser.parse_args()

    reads = fastamanipulations.fastaSubset(args.fasta, args.readIds, present=args.present)

    if not args.out:
        SeqIO.write(reads, sys.stdout, 'fasta')
    else:
        with open(args.out, 'w') as fp:
            SeqIO.write(reads, fp, 'fasta')



#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from dark import fastamanipulations


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Non-interactively generate an alignment panel',
        epilog='Given a JSON BLAST output file, a FASTA sequence file, '
        'and interesting criteria, produce an alignment panel.'
    )

    parser.add_argument(
        'json', metavar='BLAST-JSON-file', type=str, nargs='+',
        help='the JSON file of BLAST output.')

    parser.add_argument(
        '--eCutoff', type=float, default=None,
        help='Ignore hits with e-values greater (i.e., worse) than or equal '
        'to this.')

    parser.add_argument(
        '--bitCutoff', type=float, default=None,
        help='Ignore hits with bit score lower (i.e., worse) than or equal '
        'to this.')

    parser.add_argument(
        '--out', type=str, default=None,
        help='Where the output should be written to. If None, write to stdout.')

    args = parser.parse_args()

    readIds = fastamanipulations.getReadsIdsFromBlast(args.json, eCutoff=args.eCutoff, bitCutoff=args.bitCutoff)

    if not args.out:
        SeqIO.write(reads, sys.stdout, 'fasta')
    else:
        fastamanipulations.writeReadIds(readIds, args.out)


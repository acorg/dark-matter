#!/usr/bin/env python

import argparse
import sys

from dark.summarize import summarizePosition
from dark.reads import DNARead, AARead
from dark.fasta import FastaReads


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='For a fasta file with sequences, summarize what is '
                    'happening at a specific position.')

    parser.add_argument(
        '--fastaFile', required=True,
        help='The name of the FASTA file to read.')

    parser.add_argument(
        '--position', required=True, type=int,
        help='The position which should be evaluated.')

    parser.add_argument(
        '--type', required=True,
        help='The type of reads (dna or aa).')

    args = parser.parse_args()
    if args.type == 'dna':
        readClass = DNARead
    elif args.type == 'aa':
        readClass = AARead
    else:
        print >>sys.stderr, '--type must be either dna or aa.'
        sys.exit(1)

    records = FastaReads(args.fastaFile, readClass)

    result = summarizePosition(records, args.position)

    print '%d out of %d sequences were excluded, because of length.' % (
        result['includedSequences'], result['allSequences'])

    for base, count in result['countAtPosition'].items():
        print '%s: Total: %s; percentage: %s' % (
            base, count, (float(count) / float(result['includedSequences'])))

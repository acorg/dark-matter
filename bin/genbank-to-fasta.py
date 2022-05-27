#!/usr/bin/env python

import sys
import re
import argparse
from Bio import SeqIO

from dark.reads import DNARead

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Read a GenBank flat file and print FASTA for all '
                 '(or specific) records.'))

parser.add_argument(
    'genbankFile', metavar='FILE.gb', nargs='+',
    help='The GenBank file(s) to examine.')

parser.add_argument(
    '--idRegex', metavar='REGEX',
    help=('Only print records whose id (or description, if --matchDescription '
          'is used) matches this regular expression.'))

parser.add_argument(
    '--match', choices=('fullmatch', 'match', 'search'), default='search',
    help=('How to match the regular expression. See pydoc re for an '
          'explanation.'))

parser.add_argument(
    '--matchDescription', action='store_true',
    help='Match on the record description as well as the id.')

parser.add_argument(
    '--expectedCount', type=int,
    help=('The expected number of records to find. If this number is not '
          'found, print a message and exit non-zero.'))

parser.add_argument(
    '--verbose', action='store_true',
    help='Print a summary of records read and printed.')

args = parser.parse_args()

recordCount = matchCount = 0
match = getattr(re.compile(args.idRegex), args.match) if args.idRegex else None

matchDescription = args.matchDescription

for filename in args.genbankFile:
    with open(filename) as fp:
        for record in SeqIO.parse(fp, 'gb'):
            recordCount += 1
            id_ = record.id + (
                f" {record.description}" if record.description else "")
            if match is None or match(id_ if matchDescription else record.id):
                matchCount += 1
                print(DNARead(id_, str(record.seq)).toString(), end='')

if args.verbose:
    print(f'Read {recordCount} record{"s" if recordCount > 1 else ""}, '
          f'printed {matchCount}', file=sys.stderr)

if args.expectedCount is not None and matchCount != args.expectedCount:
    print(f'Expected {args.expectedCount} records but found {matchCount}. '
          f'Exiting.', file=sys.stderr)
    sys.exit(1)

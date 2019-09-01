#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Parse GenBank flat files. Exit non-zero if parsing fails.')

parser.add_argument(
    '--expectedCount', type=int,
    help=('The expected number of records in each GenBank file. If given '
          'and this number is not present in any file, exit non-zero.'))

parser.add_argument(
    'files', nargs='+',
    help='The files to parse.')

parser.add_argument(
    '--quiet', default=False, action='store_true',
    help='If given, write no output.')

args = parser.parse_args()

expectedCount = args.expectedCount
totalCount = 0

for i in args.files:
    count = 0
    try:
        records = SeqIO.parse(open(i), 'gb')
        for record in records:
            count += 1
    except ValueError:
        if not args.quiet:
            print('Could not parse %s' % i)
        sys.exit(1)
    else:
        if not args.quiet:
            print('Read %d records from %s' % (count, i))
        if expectedCount is not None and count != expectedCount:
            if not args.quiet:
                print('Expected %d records. Exiting.' % expectedCount)
            sys.exit(1)
        totalCount += count

if not args.quiet:
    print('Total records read: %d' % totalCount)

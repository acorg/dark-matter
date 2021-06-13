#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from collections import defaultdict

from dark.sam import samfile

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Print SAM/BAM file reference sequence and unmapped read '
                 'counts.'))

parser.add_argument(
    'samFile', metavar='FILENAME',
    help='The name of a SAM/BAM alignment file.')

parser.add_argument(
    '--sort', action='store_true',
    help=('Sort the output by decreasing read count (i.e., the reference with '
          'the highest number of matching reads is printed first).'))

args = parser.parse_args()


def referenceInfo():
    return {
        'duplicateCount': 0,
        'primaryCount': 0,
        'qcFailCount': 0,
        'nonDuplicateCount': 0,
        'readIds': set(),
        'secondaryCount': 0,
        'supplementaryCount': 0,
    }


referenceReads = defaultdict(referenceInfo)
mappedCount = unmappedCount = 0
readIds = set()

with samfile(args.samFile) as fp:
    for read in fp.fetch():
        readIds.add(read.query_name)
        if read.is_unmapped:
            unmappedCount += 1
        else:
            mappedCount += 1
            stats = referenceReads[read.reference_name]
            stats['readIds'].add(read.query_name)
            if read.is_secondary:
                stats['secondaryCount'] += 1
            elif read.is_supplementary:
                stats['supplementaryCount'] += 1
            else:
                stats['primaryCount'] += 1
            if read.is_duplicate:
                stats['duplicateCount'] += 1
            else:
                stats['nonDuplicateCount'] += 1
            if read.is_qcfail:
                stats['qcFailCount'] += 1

totalReads = len(readIds)

print('Found a total of %d read%s (%d mapped, %d unmapped).' %
      (totalReads, '' if totalReads == 1 else 's', mappedCount, unmappedCount))


def pct(a, b):
    return (a / b if b else 0.0) * 100.0


if args.sort:
    def key(referenceId):
        return len(referenceReads[referenceId]['readIds'])

    sortedReferenceReads = sorted(referenceReads, key=key, reverse=True)
else:
    sortedReferenceReads = sorted(referenceReads)

for referenceId in sortedReferenceReads:
    stats = referenceReads[referenceId]
    readCount = len(stats['readIds'])
    print('%s: %d/%d (%.2f%%) reads mapped to the reference.\n'
          '  Non-duplicates: %d (%.2f%%), '
          'Duplicates: %d (%.2f%%), '
          'QC fails: %d (%.2f%%)\n'
          '  Primary: %d (%.2f%%), '
          'Secondary: %d (%.2f%%), '
          'Supplementary: %d (%.2f%%)' %
          (referenceId, readCount, totalReads, pct(readCount, totalReads),
           stats['nonDuplicateCount'],
           pct(stats['nonDuplicateCount'], readCount),
           stats['duplicateCount'],
           pct(stats['duplicateCount'], readCount),
           stats['qcFailCount'],
           pct(stats['qcFailCount'], readCount),
           stats['primaryCount'],
           pct(stats['primaryCount'], readCount),
           stats['secondaryCount'],
           pct(stats['secondaryCount'], readCount),
           stats['supplementaryCount'],
           pct(stats['supplementaryCount'], readCount)))

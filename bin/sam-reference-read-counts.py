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

args = parser.parse_args()


def referenceInfo():
    return {
        'readIds': set(),
        'primaryCount': 0,
        'secondaryCount': 0,
        'supplementaryCount': 0,
    }


referenceReads = defaultdict(referenceInfo)
unmappedCount = 0
readIds = set()
mappingCount = 0

with samfile(args.samFile) as fp:
    for read in fp.fetch():
        mappingCount += 1
        readIds.add(read.query_name)
        if read.is_unmapped:
            unmappedCount += 1
        else:
            stats = referenceReads[read.reference_name]
            stats['readIds'].add(read.query_name)
            if read.is_secondary:
                stats['secondaryCount'] += 1
            elif read.is_supplementary:
                stats['supplementaryCount'] += 1
            else:
                stats['primaryCount'] += 1

totalReads = len(readIds)

print('Found a total of %d read%s, with a total of %d mapping%s and '
      '%d unmapped.' %
      (totalReads, '' if totalReads == 1 else 's',
       mappingCount, '' if mappingCount == 1 else 's', unmappedCount))

for referenceId in sorted(referenceReads):
    stats = referenceReads[referenceId]
    readCount = len(stats['readIds'])
    fraction = (readCount / totalReads) if totalReads else 0.0
    print('%s: %d (%.2f%%) reads mapped. '
          'Primary %d, secondary %d, supplementary %d' %
          (referenceId, readCount, fraction * 100.0,
           stats['primaryCount'], stats['secondaryCount'],
           stats['supplementaryCount']))

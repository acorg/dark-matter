#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from collections import defaultdict

from dark.sam import samfile
from dark.utils import pct


def referenceInfo():
    return {
        'duplicate': set(),
        'primary': set(),
        'qcFail': set(),
        'nonDuplicate': set(),
        'readIds': set(),
        'secondary': set(),
        'supplementary': set(),
    }


def main(args):
    referenceReads = defaultdict(referenceInfo)
    mapped = set()
    unmapped = set()
    readIds = set()

    with samfile(args.samFile) as fp:
        for read in fp.fetch():
            id_ = read.query_name
            readIds.add(id_)
            if read.is_unmapped:
                unmapped.add(id_)
            else:
                mapped.add(id_)
                stats = referenceReads[read.reference_name]
                stats['readIds'].add(id_)

                if read.is_secondary:
                    stats['secondary'].add(id_)
                elif read.is_supplementary:
                    stats['supplementary'].add(id_)
                else:
                    stats['primary'].add(id_)

                if read.is_duplicate:
                    stats['duplicate'].add(id_)
                else:
                    stats['nonDuplicate'].add(id_)

                if read.is_qcfail:
                    stats['qcFail'].add(id_)

    totalReads = len(readIds)

    print('Found a total of %d read%s (%d mapped, %d unmapped).' %
          (totalReads, '' if totalReads == 1 else 's',
           len(mapped), len(unmapped)))

    def key(referenceId):
        return len(referenceReads[referenceId]['readIds'])

    sortedReferenceReads = sorted(referenceReads, key=key, reverse=True)
    topReference = sortedReferenceReads[0]

    if not args.sort:
        # Just sort the references by name, not by number of reads.
        sortedReferenceReads = sorted(referenceReads)

    cumulativeReadIds = set()

    for referenceId in sortedReferenceReads:
        stats = referenceReads[referenceId]
        readCount = len(stats['readIds'])
        newReadCount = len(stats['readIds'] - cumulativeReadIds)
        cumulativeReadIds.update(stats['readIds'])
        print('\n%s:\n'
              '  Overall reads mapped to the reference: %s\n'
              '  Non-duplicates: %s, Duplicates: %s, QC fails: %s\n'
              '  Primary: %s, Secondary: %s, Supplementary: %s\n'
              '  Reads not matching any reference above: %s' %
              (referenceId,
               pct(readCount, totalReads),
               pct(len(stats['nonDuplicate']), readCount),
               pct(len(stats['duplicate']), readCount),
               pct(len(stats['qcFail']), readCount),
               pct(len(stats['primary']), readCount),
               pct(len(stats['secondary']), readCount),
               pct(len(stats['supplementary']), readCount),
               pct(newReadCount, totalReads)))

    # Write out the (sorted) read ids of the reference with the most reads.
    if args.topReferenceIdsFile:
        with open(args.topReferenceIdsFile, 'w') as fp:
            print('\n'.join(sorted(referenceReads[topReference]['readIds'])),
                  file=fp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Print SAM/BAM file reference sequence and unmapped read '
                     'counts.'))

    parser.add_argument(
        'samFile', metavar='FILENAME',
        help='The name of a SAM/BAM alignment file.')

    parser.add_argument(
        '--sort', action='store_true',
        help=('Sort the output by decreasing read count (i.e., the reference '
              'with the highest number of matching reads is printed first).'))

    parser.add_argument(
        '--topReferenceIdsFile', metavar='FILE',
        help=('The file to write the (sorted) ids of the reads for the '
              'reference with the highest number of matching reads.'))

    main(parser.parse_args())

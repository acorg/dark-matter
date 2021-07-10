#!/usr/bin/env python

from __future__ import print_function, division

import sys
import argparse
from collections import defaultdict

from dark.sam import samfile, SAMFilter
from dark.utils import pct


def referenceInfo():
    """
    Hold information about the reads that match a reference.

    @return: A C{dict}, as below.
    """
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
    """
    Print SAM/BAM file reference read counts.

    @param args: An argparse namespace with information about parsed
        command-line options.
    """
    if args.topReferenceIdsFile and args.sortBy != 'count':
        print('--topReferenceIdsFile only makes sense when using --sortBy '
              'count', file=sys.stderr)
        sys.exit(1)

    referenceReads = defaultdict(referenceInfo)
    mapped = set()
    unmapped = set()
    readIds = set()

    referenceLengths = SAMFilter(args.samFile).referenceLengths()

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

    print('Found a total of %d read%s (%d mapped, %d unmapped) mapping '
          'against %d of %d reference%s.' %
          (totalReads, '' if totalReads == 1 else 's',
           len(mapped), len(unmapped),
           len(referenceReads), len(referenceLengths),
           '' if len(referenceLengths) == 1 else 's'))

    if args.sortBy == 'count':

        def key(referenceId):
            return len(referenceReads[referenceId]['readIds'])

        sortedReferenceReads = sorted(referenceReads, key=key, reverse=True)
        topReference = sortedReferenceReads[0]
    else:
        # Sort the references by name
        sortedReferenceReads = sorted(referenceReads)

    cumulativeReadIds = set()

    for count, referenceId in enumerate(sortedReferenceReads, start=1):
        stats = referenceReads[referenceId]
        readCount = len(stats['readIds'])
        if readCount == 0 and args.excludeZeroes:
            continue
        newReadCount = len(stats['readIds'] - cumulativeReadIds)
        if newReadCount == 0 and args.excludeIfNoAdditional:
            continue
        cumulativeReadIds.update(stats['readIds'])
        print('\nReference %d: %s (%d nt):\n'
              '  Overall reads mapped to the reference: %s\n'
              '  Non-duplicates: %s, Duplicates: %s, QC fails: %s\n'
              '  Primary: %s, Secondary: %s, Supplementary: %s\n'
              '  Reads not matching any reference above: %s\n'
              '  Previously unmatched reads for this reference: %s' %
              (count, referenceId, referenceLengths[referenceId],
               pct(readCount, totalReads),
               pct(len(stats['nonDuplicate']), readCount),
               pct(len(stats['duplicate']), readCount),
               pct(len(stats['qcFail']), readCount),
               pct(len(stats['primary']), readCount),
               pct(len(stats['secondary']), readCount),
               pct(len(stats['supplementary']), readCount),
               pct(newReadCount, totalReads),
               pct(newReadCount, readCount)))

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
        '--sortBy', choices=('count', 'name'), default='count',
        help=('Use "count" to sort the output by decreasing read count (i.e., '
              'the reference with the highest number of matching reads is '
              'printed first), or "name" to sort by reference name.'))

    parser.add_argument(
        '--excludeZeroes', action='store_true',
        help='Do not print references unmatched by any read.')

    parser.add_argument(
        '--excludeIfNoAdditional', action='store_true',
        help=('Do not print references that are mapped by no additional reads '
              '(i.e., reads other than those matching already-printed '
              'references). This is useful when --sortBy count is used, '
              'making it possible to not print references only matched by '
              'reads matching references that have already been printed.'))

    parser.add_argument(
        '--topReferenceIdsFile', metavar='FILE',
        help=('The file to write the (sorted) ids of the reads for the '
              'reference with the highest number of matching reads. Only '
              'valid when --sortBy count is used.'))

    main(parser.parse_args())

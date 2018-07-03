#!/usr/bin/env python

from __future__ import print_function, division

import sys
from pysam import AlignmentFile

from dark.filter import (
    addFASTAFilteringCommandLineOptions, parseFASTAFilteringCommandLineOptions)
from dark.reads import Read, Reads
from dark.sam import samfile


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Given a SAM/BAM file and a set of filtering criteria '
                     'write filtered SAM/BAM to stdout.'))

    parser.add_argument(
        'samfile',
        help='The SAM/BAM file to filter.')

    parser.add_argument(
        '--quiet', action='store_true', default=False,
        help='If True, do not print the final summary.')

    parser.add_argument(
        '--bam', action='store_const', const='b', default='',
        help='If given, write (gzip compressed) BAM output.')

    parser.add_argument(
        '--checkResultCount', type=int,
        help=('The number of alignments expected in the output. If this '
              'number is not seen, the script exits with status 1 (and an '
              'error message is printed unless --quiet is used).'))

    parser.add_argument(
         '--dropUnmapped', default=False, action='store_true',
         help='If given, unmapped matches will not be output.')

    parser.add_argument(
         '--dropSecondary', default=False, action='store_true',
         help='If given, secondary matches will not be output.')

    parser.add_argument(
         '--dropSupplementary', default=False, action='store_true',
         help='If given, supplementary matches will not be output.')

    parser.add_argument(
         '--dropDuplicates', default=False, action='store_true',
         help=('If given, matches flagged as optical or PCR duplicates will '
               'not be output.'))

    parser.add_argument(
         '--keepQCFailures', default=False, action='store_true',
         help=('If given, reads that are considered quality control failures '
               'will be included in the output.'))

    parser.add_argument(
        '--referenceWhitelist', metavar='NAME', action='append',
        help=('A reference sequence id whose alignments should be output. '
              'If omitted, alignments against all references will be output. '
              'May be repeated.'))

    parser.add_argument(
        '--referenceBlacklist', metavar='NAME', action='append',
        help=('A reference sequence id whose alignments should not be output. '
              'If omitted, alignments against all references will be output. '
              'May be repeated.'))

    addFASTAFilteringCommandLineOptions(parser, filteringSAM=True)

    args = parser.parse_args()

    referenceWhitelist = (set(args.referenceWhitelist)
                          if args.referenceWhitelist else None)
    referenceBlacklist = (set(args.referenceBlacklist)
                          if args.referenceBlacklist else None)

    if referenceWhitelist and referenceBlacklist and (
            referenceWhitelist & referenceBlacklist):
        raise ValueError(
            'The reference whitelist and blacklist sets are not '
            'disjoint. The intersection contains %s.' %
            ', '.join(sorted(referenceWhitelist & referenceBlacklist)))

    reads = parseFASTAFilteringCommandLineOptions(args, Reads(),
                                                  filteringSAM=True)
    filterRead = reads.filterRead
    dropUnmapped = args.dropUnmapped
    dropSecondary = args.dropSecondary
    dropSupplementary = args.dropSupplementary
    dropDuplicates = args.dropDuplicates
    keepQCFailures = args.keepQCFailures

    kept = total = 0
    with samfile(args.samfile) as samAlignment:
        out = AlignmentFile(sys.stdout, mode='w' + args.bam,
                            template=samAlignment)
        save = out.write
        for alignment in samAlignment.fetch():
            total += 1

            if (filterRead(Read(alignment.query_name,
                                alignment.query_sequence,
                                alignment.qual))
                    and not (
                        (alignment.is_unmapped and dropUnmapped) or
                        (alignment.is_secondary and dropSecondary) or
                        (alignment.is_supplementary and dropSupplementary) or
                        (alignment.is_duplicate and dropDuplicates) or
                        (alignment.is_qcfail and not keepQCFailures) or
                        (referenceWhitelist is not None and
                         alignment.reference_name not in referenceWhitelist) or
                        (referenceBlacklist is not None and
                         alignment.reference_name in referenceBlacklist))):
                kept += 1
                save(alignment)
        out.close()

    if not args.quiet:
        print('Read %d alignment%s, kept %d (%.2f%%).' %
              (total, '' if total == 1 else 's', kept,
               0.0 if total == 0 else kept / total * 100.0), file=sys.stderr)

    if args.checkResultCount is not None:
        if kept != args.checkResultCount:
            if not args.quiet:
                print('Did not write the expected %d alignment%s (wrote %d).' %
                      (args.checkResultCount,
                       '' if args.checkResultCount == 1 else 's', kept),
                      file=sys.stderr)
            sys.exit(1)

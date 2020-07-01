#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from collections import Counter
from numpy import std

from dark.filter import (
    addFASTAFilteringCommandLineOptions, parseFASTAFilteringCommandLineOptions)
from dark.reads import Reads
from dark.sam import samfile, SAMFilter, samReferences, UnknownReference
from dark.utils import baseCountsToStr, pct


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Print SAM/BAM file coverage statistics by offset. '
                 'Output lines show the offset.'))

addFASTAFilteringCommandLineOptions(parser)
SAMFilter.addFilteringOptions(parser, samfileIsPositional=True)

parser.add_argument(
    '--noOffsets', default=False, action='store_true',
    help='Do not print per-offset details of base counts.')

parser.add_argument(
    '--noStats', default=False, action='store_true',
    help='Do not print final average and standard deviation statistics.')

parser.add_argument(
    '--noFilter', default=False, action='store_true',
    help=('Do not use our SAM filtering. Note that if you give this option, '
          'any filtering option (other than --referenceId) you also specify '
          'that is provided by the SAMFilter.addFilteringOptions will be '
          'silently ignored!'))

args = parser.parse_args()

if args.noOffsets and args.noStats:
    print('You have used both --noOffsets and --noStats, so there is no '
          'output!', file=sys.stderr)
    sys.exit(1)


# We don't have a file of reads, we just want a read filter that we can use
# to filter the SAM file query sequences and to get reference lengths from.
reads = parseFASTAFilteringCommandLineOptions(args, Reads())
samFilter = SAMFilter.parseFilteringOptions(args, reads.filterRead)

printOffsets = not args.noOffsets
printStats = not args.noStats

if samFilter.referenceIds and len(samFilter.referenceIds) > 1:
    print('Only one reference id can be given. To calculate coverage for more '
          'than one reference, run this script multiple times.',
          file=sys.stderr)
    sys.exit(1)

try:
    referenceLengths = samFilter.referenceLengths()
except UnknownReference:
    referenceId = samFilter.referenceIds.pop()
    referenceIds = samReferences(args.samfile)
    print('Reference %r does not appear in SAM file %s. Known '
          'references are: %s.' % (
              referenceId, args.samfile, ', '.join(sorted(referenceIds))),
          file=sys.stderr)
    sys.exit(1)


if args.noFilter:
    # Do not do our custom SAM filtering.
    def filterRead(read):
        return not (read.is_del or read.is_refskip)
else:
    def filterRead(read):
        return (not (read.is_del or read.is_refskip) and
                samFilter.filterAlignment(read.alignment))


if printStats:
    counts = []

with samfile(args.samfile) as sam:

    if samFilter.referenceIds:
        # No need to check if the given reference id is in referenceLengths
        # because the samFilter.referenceLengths call above catches that.
        referenceId = samFilter.referenceIds.pop()
    else:
        if len(referenceLengths) == 1:
            referenceId = list(referenceLengths)[0]
        else:
            print('SAM file %r contains %d references (%s). Only one '
                  'reference id can be analyzed at a time. Please use '
                  '--referenceId to specify the one you want examined.' % (
                      args.samfile, len(referenceLengths),
                      ', '.join(sorted(referenceLengths))), file=sys.stderr)
            sys.exit(1)

    for column in sam.pileup(reference=referenceId):
        bases = Counter()
        for read in column.pileups:
            if filterRead(read):
                base = read.alignment.query_sequence[read.query_position]
                bases[base] += 1

        baseCount = sum(bases.values())

        if printStats:
            counts.append(baseCount)

        if printOffsets:
            print('%d: %d %s' % (column.reference_pos + 1, baseCount,
                                 baseCountsToStr(bases)))

if printStats:
    referenceLength = referenceLengths[referenceId]
    print('Reference id: %s' % referenceId)
    print('Reference length: %d' % referenceLength)
    print('Bases covered: %s' % pct(len(counts), referenceLength))
    print('Min coverage depth: %d' % (
        0 if len(counts) < referenceLength else min(counts)))
    if counts:
        # Don't use Python3 default= option on max. Trying to keep Python 2
        # compatibility.
        print('Max coverage depth: %d' % max(counts))
    else:
        print('Max coverage depth: 0')
    print('Mean coverage depth: %.3f' % (sum(counts) / referenceLength))
    if counts:
        print('Coverage depth s.d.: %.3f' % std(counts))

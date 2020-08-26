#!/usr/bin/env python

from __future__ import print_function

import argparse
from collections import defaultdict, Counter

from dark.filter import (
    addFASTAFilteringCommandLineOptions, parseFASTAFilteringCommandLineOptions)
from dark.reads import Reads
from dark.sam import samfile, SAMFilter

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Print SAM/BAM file coverage statistics.')

parser.add_argument(
    '--noFilter', default=False, action='store_true',
    help=('Do not use our SAM filtering. Note that if you give this option, '
          'any filtering option (other than --referenceId) you also specify '
          'that is provided by the SAMFilter.addFilteringOptions will be '
          'silently ignored!'))

addFASTAFilteringCommandLineOptions(parser)
SAMFilter.addFilteringOptions(parser, samfileIsPositional=True)

args = parser.parse_args()

if args.noFilter:
    # Do not do our custom SAM filtering.
    def filterRead(read):
        return True
else:
    def filterRead(read):
        return (not (read.is_del or read.is_refskip) and
                samFilter.filterAlignment(read.alignment))

# We don't have a file of reads, we just want a read filter that we can use
# to filter the SAM file query sequences and to get reference lengths from.
reads = parseFASTAFilteringCommandLineOptions(args, Reads())
samFilter = SAMFilter.parseFilteringOptions(args, reads.filterRead)

coveredOffsets = defaultdict(Counter)
coveringReads = defaultdict(set)

with samfile(args.samfile) as sam:
    for column in sam.pileup():
        referenceOffset = column.reference_pos
        for read in column.pileups:
            if filterRead(read):
                referenceId = read.alignment.reference_name
                coveredOffsets[referenceId][referenceOffset] += 1
                coveringReads[referenceId].add(read.alignment.query_name)

referenceLengths = samFilter.referenceLengths()

for referenceId in sorted(referenceLengths):
    offsetsCovered = len(coveredOffsets[referenceId])
    referenceLength = referenceLengths[referenceId]
    print('%s: length %d, covering reads %d, covered sites %d (%.4f%%), '
          'mean coverage depth %.4f (min: %d, max: %d)' %
          (referenceId, referenceLength, len(coveringReads[referenceId]),
           offsetsCovered, offsetsCovered / referenceLength * 100.0,
           sum(coveredOffsets[referenceId].values()) / referenceLength,
           min(coveredOffsets[referenceId].values(), default=0),
           max(coveredOffsets[referenceId].values(), default=0)))

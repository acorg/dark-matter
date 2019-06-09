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

addFASTAFilteringCommandLineOptions(parser)
SAMFilter.addFilteringOptions(parser, samfileIsPositionalArg=True)

args = parser.parse_args()

# We don't have a file of reads, we just want a read filter that we
# can use to filter the SAM file query sequences.
reads = parseFASTAFilteringCommandLineOptions(args, Reads())
samFilter = SAMFilter.parseFilteringOptions(args, reads.filterRead)

coveredOffsets = defaultdict(Counter)
coveringReads = defaultdict(set)

with samfile(args.samfile) as sam:
    for column in sam.pileup():
        referenceOffset = column.reference_pos
        for read in column.pileups:
            if samFilter.filterAlignment(read.alignment):
                referenceId = read.alignment.reference_name
                coveredOffsets[referenceId][referenceOffset] += 1
                coveringReads[referenceId].add(read.alignment.query_name)

referenceLengths = samFilter.referenceLengths()

for referenceId in sorted(referenceLengths):
    offsetsCovered = len(coveredOffsets[referenceId])
    referenceLength = referenceLengths[referenceId]
    print('%s: length %d, covering reads %d, covered sites %d (%.4f%%), '
          'coverage depth %.4f' %
          (referenceId, referenceLength, len(coveringReads[referenceId]),
           offsetsCovered, offsetsCovered / referenceLength * 100.0,
           sum(coveredOffsets[referenceId].values()) / referenceLength))

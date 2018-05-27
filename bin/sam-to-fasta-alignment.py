#!/usr/bin/env python

"""
Convert SAM/BAM to an alignment in FASTA or Phylip format.
"""

from __future__ import division, print_function

import sys
import argparse
from dark.reads import Read
from pysam import (
    AlignmentFile, CMATCH, CINS, CDEL, CREF_SKIP, CSOFT_CLIP, CHARD_CLIP, CPAD,
    CEQUAL, CDIFF)

# See https://samtools.github.io/hts-specs/SAMv1.pdf for the True/False values
# in the following and http://pysam.readthedocs.io/en/latest/api.html for the
# pysam documentation.

CONSUMES_QUERY = {
    int(CMATCH): True,
    int(CINS): True,
    int(CDEL): False,
    int(CREF_SKIP): False,
    int(CSOFT_CLIP): True,
    int(CHARD_CLIP): False,
    int(CPAD): False,
    int(CEQUAL): True,
    int(CDIFF): True,
}

CONSUMES_REFERENCE = {
    int(CMATCH): True,
    int(CINS): False,
    int(CDEL): True,
    int(CREF_SKIP): True,
    int(CSOFT_CLIP): False,
    int(CHARD_CLIP): False,
    int(CPAD): False,
    int(CEQUAL): True,
    int(CDIFF): True,
}

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Convert a SAM/BAM file to aligned FASTA.'))

parser.add_argument(
    'samFile', help='The SAM/BAM file to convert.')

parser.add_argument(
     '--minLength', type=int, default=-1,
     help='Do not output anything for queries shorter than this.')


args = parser.parse_args()
samfile = AlignmentFile(args.samFile)

if samfile.nreferences != 1:
    print('The SAM/BAM file must have only 1 reference sequence (yours has '
          ' %d).' % samfile.nreferences, file=sys.stderr)
    sys.exit(1)

referenceLength = samfile.lengths[0]
minLength = args.minLength

for read in samfile.fetch():
    if read.is_unmapped:
        continue
    query = read.query_sequence
    queryLength = len(query)
    if queryLength < minLength:
        continue
    index = 0
    alignedSequence = ['-'] * read.reference_start
    for operation, length in read.cigartuples:
        operation = int(operation)
        if CONSUMES_QUERY[operation]:
            alignedSequence.append(query[index:index + length])
            index += length
        elif CONSUMES_REFERENCE[operation]:
            alignedSequence.append('-' * length)

    # Sanity check that we consumed the entire query.
    assert index == queryLength

    alignedSequence = ''.join(alignedSequence)
    alignedSequence += '-' * (referenceLength - len(alignedSequence))

    print(Read(read.query_name, alignedSequence).toString('fasta'), end='')

samfile.close()

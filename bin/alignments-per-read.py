#!/usr/bin/env python

"""
Read in BLAST records and print a count of how many sequences
each read matches with, followed by the name of the read.
"""

from __future__ import print_function

import sys

from dark.blast.alignments import BlastReadsAlignments
from dark.fasta import FastaReads

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: %s reads.fasta BLAST-files.json...' %
              sys.argv[0], file=sys.stderr)
    else:
        readsFile = sys.argv[1]
        jsonFiles = sys.argv[2:]
        reads = FastaReads(readsFile)
        readsAlignments = BlastReadsAlignments(readsFile, jsonFiles)
        for readAlignments in readsAlignments:
            print(len(readAlignments.alignments), readAlignments.read.id)

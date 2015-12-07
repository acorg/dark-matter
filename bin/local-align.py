#!/usr/bin/env python

"""
Aligns each sequence in fileHandle1 to each sequence in fileHandle2.
"""

from __future__ import print_function

from dark.local_align import LocalAlignment
import sys
from Bio import SeqIO

if len(sys.argv) != 3 or len(sys.argv) != 8:
    print('Usage: %s filename1, filename2, (matchScore, '
          'mismatchScore, gapOpenScore, gapExtendScore, '
          'gapExtendDecay)'
          % sys.argv[0], file=sys.stderr)
    sys.exit(1)
else:
    fileHandle1 = sys.argv[1]
    fileHandle2 = sys.argv[2]
    if len(sys.argv) == 3:
        matchScore = 1
        mismatchScore = -1
        gapOpenScore = -1
        gapExtendScore = -1
        gapExtendDecay = 0.0
    else:
        matchScore = int(sys.argv[3])
        mismatchScore = int(sys.argv[4])
        gapOpenScore = int(sys.argv[5])
        gapExtendScore = int(sys.argv[6])
        gapExtendDecay = float(sys.argv[7])
    for record in SeqIO.parse(fileHandle1, 'fasta'):
        input1 = record
        for record in SeqIO.parse(fileHandle2, 'fasta'):
            input2 = record
            alignment = LocalAlignment(input1, input2, match=matchScore,
                                       mismatch=mismatchScore,
                                       gap=gapOpenScore,
                                       gapExtend=gapExtendScore,
                                       gapDecay=gapExtendDecay)
            print(alignment.createAlignment())

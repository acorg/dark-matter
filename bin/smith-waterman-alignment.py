#!/usr/bin/env python

"""
Aligns each sequence in fileHandle1 to each sequence in fileHandle2.
Returns list of three strings for each alignment.
"""

from dark.smith_waterman import smith_waterman
import sys
from Bio import SeqIO

if len(sys.argv) < 3 or len(sys.argv) > 6:
    print >>sys.stderr, "Usage: %s filename1, filename2, (matchScore, mismatchScore, gapScore)" % sys.argv[0]
    sys.exit(1)

else:
    fileHandle1 = sys.argv[1]
    fileHandle2 = sys.argv[2]
    if len(sys.argv) == 3:
        matchScore = 1
        mismatchScore = -1
        gapScore = -1
    else:
        matchScore = int(sys.argv[3])
        mismatchScore = int(sys.argv[4])
        gapScore = int(sys.argv[5])
    for record in SeqIO.parse(fileHandle1, "fasta"):
        input1 = record.seq
        for record in SeqIO.parse(fileHandle2, "fasta"):
            input2 = record.seq
            result = smith_waterman(input1, input2, match=matchScore, mismatch=mismatchScore, gap=gapScore)
            print "\n".join(result)

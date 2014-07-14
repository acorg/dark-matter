#!/usr/bin/env python

"""
Aligns each sequence in fileHandle1 to each sequence in fileHandle2. 
Returns list of three strings for each alignment.
"""

from dark.smith_waterman import smith_waterman
import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    print >>sys.stderr, "Usage: %s filename1, filename2" % sys.argv[0]
    sys.exit(1)

else:
    fileHandle1 = sys.argv[1]
    fileHandle2 = sys.argv[2]
    for record in SeqIO.parse(fileHandle1, "fasta"):
        input1 = record.seq.upper()
        for record in SeqIO.parse(fileHandle2, "fasta"):
            input2 = record.seq.upper()
            result = smith_waterman(input1, input2)
            print "\n".join(result)

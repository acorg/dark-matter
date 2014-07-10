#!/usr/bin/env python

from dark.smith_waterman import smith_waterman
import sys
from Bio import SeqIO

if len(sys.argv) > 3:
    print >>sys.stderr, "Usage: %s" % sys.argv[0]
    sys.exit(1)

else:
	file_handle1 = sys.argv[1]
	file_handle2 = sys.argv[2]
	for record in SeqIO.parse(file_handle1, "fasta"):
		input1 = record.seq.upper()
		for record in SeqIO.parse(file_handle2, "fasta"):
			input2 = record.seq.upper()
			result = smith_waterman(input1, input2)
			print "\n".join(result)


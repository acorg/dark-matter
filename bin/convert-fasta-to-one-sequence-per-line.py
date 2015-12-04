#!/usr/bin/env python

"""
Read FASTA records from stdin or a given file name, write the same
records to stdout but with sequences taking up only one line. This is to
reformat FASTA files into something that can be more easily processed by
UNIX tools (like split).
"""

from __future__ import print_function

from Bio import SeqIO
import sys

if len(sys.argv) == 1:
    fp = sys.stdin
else:
    fp = open(sys.argv[1])

for seq in SeqIO.parse(fp, 'fasta'):
    print('>%s\n%s' % (seq.description, str(seq.seq)))

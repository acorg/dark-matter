#!/usr/bin/env python

from Bio import SeqIO
import sys

if len(sys.argv) < 2:
    print >>sys.stderr, 'Usage: %s file.sff > file.fastq' % sys.argv[0]
    sys.exit(1)
else:
    SeqIO.convert(sys.argv[1], 'sff', sys.stdout, 'fastq')

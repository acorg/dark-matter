#!/usr/bin/env python

import sys

from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: %s file.sff > file.fastq" % sys.argv[0], file=sys.stderr)
    sys.exit(1)
else:
    SeqIO.convert(sys.argv[1], "sff", sys.stdout, "fastq")

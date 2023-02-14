#!/usr/bin/env python

"""
Read DNA FASTA from stdin and print FASTA to stdout, with the sequences
randomized. The read ids and lengths are preserved.

Note: This produces DNA sequences. If you have AA reads and you need this
      functionality, we can add it.
"""

import sys
from numpy.random import choice

from dark.fasta import FastaReads
from dark.reads import DNARead


if __name__ == "__main__":
    for read in FastaReads(sys.stdin):
        seq = "".join(choice(["A", "C", "G", "T"], len(read), replace=True))
        print(DNARead(read.id, seq).toString("fasta"))

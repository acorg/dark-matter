#!/usr/bin/env python

from __future__ import print_function

import sys
from Bio import SeqIO
from dark import fasta

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Usage: %s 1.fasta, 2.fasta, ... > seq.fasta' % (
            sys.argv[0]), file=sys.stderr)
        sys.exit(1)
    else:
        reads = fasta.fastaSubtract(map(open, sys.argv[1:]))
        SeqIO.write(reads, sys.stdout, 'fasta')

#!/usr/bin/env python

import sys
from Bio import SeqIO
from dark import fasta

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print >>sys.stderr, 'Usage: %s 1.fasta, 2.fasta, ... > seq.fasta' % (
            sys.argv[0])
        sys.exit(1)
    else:
        reads = fasta.fastaSubtract(sys.argv[1:])

        SeqIO.write(reads, sys.stdout, 'fasta')
        print >>sys.stderr, 'Wrote %d sequences.' % len(reads)

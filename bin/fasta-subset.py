#!/usr/bin/env python

"""
Get a set of FASTA sequence identifiers from sys.argv, read FASTA from
stdin, and print FASTA to stdout for the given sequences.
"""

import sys
from Bio import SeqIO

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print >>sys.stderr, 'Usage: %s seqId-1 seqId-2 ... < sequences.fasta'
        sys.exit(1)

    wanted = set(sys.argv[1:])
    found = []
    for seq in SeqIO.parse(sys.stdin, 'fasta'):
        if seq.description in wanted:
            wanted.remove(seq.description)
            found.append(seq)

    if found:
        SeqIO.write(found, sys.stdout, 'fasta')

    print >>sys.stderr, 'Found %d sequences.' % len(found)
    if wanted:
        print >>sys.stderr, 'WARNING: %d sequence%s not found: %s' % (
            len(wanted), '' if len(wanted) == 1 else 's were',
            ' '.join(wanted))

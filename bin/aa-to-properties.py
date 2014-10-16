#!/usr/bin/env python

"""
Read AA FASTA from stdin and print properties FASTA to stdout.
"""

import sys

from dark.fasta import FastaReads


if __name__ == '__main__':
    if len(sys.argv) > 1:
        print >>sys.stderr, ('Usage: %s < aaFile > propertiesFile.fasta' %
                             sys.argv[0])
        sys.exit(1)

    reads = FastaReads(sys.stdin)
    write = sys.stdout.write

    for read in reads:
        for aa in read.aaToProperties():
            write(aa.toString('fasta'))

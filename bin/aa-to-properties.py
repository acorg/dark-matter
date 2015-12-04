#!/usr/bin/env python

"""
Read AA FASTA from stdin and print properties FASTA to stdout.
"""

from __future__ import print_function

import sys

from dark.fasta import FastaReads


if __name__ == '__main__':
    if len(sys.argv) > 1:
        print('Usage: %s < aaFile > propertiesFile.fasta' %
              sys.argv[0], file=sys.stderr)
        sys.exit(1)

    reads = FastaReads(sys.stdin)
    write = sys.stdout.write

    for read in reads:
        for aa in read.aaToProperties():
            write(aa.toString('fasta'))

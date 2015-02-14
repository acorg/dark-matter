#!/usr/bin/env python

import sys
from os.path import basename

from dark.fastq import FastqReads

# Print a usage message if any arguments were given on the command
# line. This is to remind people who want to provide filenames that they
# must use '<' to supply our stdin and '>' to store our output.

if len(sys.argv) > 1:
    print >>sys.stderr, (
        'Usage: %s <input.fastq [>output.fasta]' % basename(sys.argv[0]))
    sys.exit(1)
else:
    write = sys.stdout.write
    for read in FastqReads(sys.stdin):
        write(read.toString('fasta'))

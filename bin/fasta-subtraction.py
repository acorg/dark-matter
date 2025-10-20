#!/usr/bin/env python

import sys

from dark import fasta

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(
            "Usage: %s 1.fasta, 2.fasta, ... > seq.fasta" % (sys.argv[0]),
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        for read in fasta.fastaSubtract(sys.argv[1:]):
            print(read.toString(), end="")

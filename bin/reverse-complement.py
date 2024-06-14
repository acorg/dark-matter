#!/usr/bin/env python

import sys

from dark.reads import DNARead


if len(sys.argv) < 2:
    print("Usage: reverse-complement.py SEQ1 [SEQ2, SEQ3, ...]", file=sys.stderr)
    sys.exit(1)
else:
    for seq in sys.argv[1:]:
        print(DNARead("id", seq).reverseComplement().sequence)

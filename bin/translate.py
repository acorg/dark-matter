#!/usr/bin/env python

import sys
from Bio.Seq import Seq

if len(sys.argv) == 0:
    print(f"Usage: {sys.argv[0]} DNA-sequence", file=sys.stderr)
    sys.exit(1)

for dna in sys.argv[1:]:
    print(Seq(dna).translate())

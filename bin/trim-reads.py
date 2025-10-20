#!/usr/bin/env python

import sys

from dark.analyze_reads import getPrefixAndSuffix, trimReads
from dark.reads import Reads

if len(sys.argv) != 2:
    sys.exit(f"{sys.argv[0]} takes one FASTA filename argument")

else:
    filename = sys.argv[1]
    with open(filename) as fp:
        prefix, suffix = getPrefixAndSuffix(fp)

    print(f"Prefix length {prefix}, suffix length {suffix}", file=sys.stderr)

    with open(filename) as fp:
        Reads(trimReads(prefix, suffix, fp)).save(sys.stdout)

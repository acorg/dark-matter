#!/usr/bin/env python

import sys

from Bio import SeqIO

from dark.analyze_reads import getPrefixAndSuffix, trimReads

if len(sys.argv) != 2:
    print("getPrefixAndSuffix() takes exactly 1 argument", file=sys.stderr)
    sys.exit(1)

else:
    filename = sys.argv[1]
    prefix, suffix = getPrefixAndSuffix(filename)

    print("Prefix length %d, suffix length %d" % (prefix, suffix), file=sys.stderr)

    SeqIO.write(trimReads(prefix, suffix, filename), sys.stdout, "fasta")

#!/usr/bin/env python

from Bio import SeqIO
from dark.analyze_reads import getPrefixAndSuffix
from dark.analyze_reads import trimReads
import sys


if len(sys.argv) != 2:
    print >> sys.stderr, "getPrefixAndSuffix() takes exactly 1 argument"
    sys.exit(1)

else:
    filename = sys.argv[1]
    prefix, suffix = getPrefixAndSuffix(filename)

    print >>sys.stderr, "Prefix length %d, suffix length %d" % (prefix, suffix)

    SeqIO.write(trimReads(prefix, suffix, filename), sys.stdout, 'fasta')

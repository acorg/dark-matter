from dark.analyze_reads import getReads, _longestPrefixOfTwoSeqs
from Bio import SeqIO
import sys


if len(sys.argv) > 2:
    print >> stderr, "ERROR, takes at least two arguments."
    sys.exit(1)

else:
    caller = open(sys.argv[1], "rU")
    getReads(caller)



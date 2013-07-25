from dark.analyze_reads import getPrefixAndSuffix
from dark.analyze_reads import trimReads
import sys


if len(sys.argv) != 2:
    print >> sys.stderr, "getPrefixAndSuffix() takes exactly 1 argument"
    sys.exit(1)

else:
    filename = sys.argv[1]
    prefix, suffix = getPrefixAndSuffix(filename)

    print "prefix %d, suffix %d" % (prefix, suffix)

    trimReads(prefix, suffix, filename)

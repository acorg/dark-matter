from dark.analyze_reads import getPrefixAndSuffix, trimReads
import sys


if len(sys.argv) > 2:
    print >> sys.stderr, "ERROR, takes at least two arguments."
    sys.exit(1)

else:
    filename = sys.argv[1]
    result = getPrefixAndSuffix(filename)
    prefix = result[0]
    suffix = result[1]

    print "prefix", prefix
    print "suffix", suffix

    trimReads(prefix, suffix, filename)

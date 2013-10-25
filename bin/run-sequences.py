#!/bin/python

from dark.sequences import summarizeRecordsBySequence
from time import time
import sys

if len(sys.argv) > 4:
    print >>sys.stderr, "Usage: %s file.fasta" % sys.argv[0]
    sys.exit(1)

else:
    filename = sys.argv[1]
    start = time()
    result = summarizeRecordsBySequence(filename)
    elapsed = time() - start

    resultdict = result[0]
    count = result[1]

    print 'Read %d BLAST records in %.3f secs (%.0f records/sec)' % (
        count, elapsed, float(count) / float(elapsed))


#!/usr/bin/env python

"""
Read simplified JSON BLAST records and report the elapsed time.
"""

from dark.conversion import readJSONRecords
from time import time
import sys


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Usage: %s file.json' % sys.argv[0]
        sys.exit(1)
    else:
        start = time()
        for count, record in enumerate(readJSONRecords(sys.argv[1])):
            pass
        stop = time()
        print 'Read %d JSON BLAST records in %.3f seconds' % (
            count, stop - start)

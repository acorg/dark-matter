#!/usr/bin/env python

"""
Read in BLAST records and print a count of how many sequences
each read matches with, followed by the name of the read.
"""

from dark.utils import readBlastRecords

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Usage: %s BLAST-file.{json,xml}' % sys.argv[0]
    else:
        results = readBlastRecords(sys.argv[1])
        for record in results:
            print len(record.alignments), record.query


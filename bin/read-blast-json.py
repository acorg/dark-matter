#!/usr/bin/env python

"""
Read simplified JSON BLAST records and report the elapsed time.
"""

from __future__ import print_function

from dark.conversion import JSONRecordsReader
from time import time
import sys


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s file.json' % sys.argv[0], file=sys.stderr)
        sys.exit(1)
    else:
        start = time()
        jsonReader = JSONRecordsReader(sys.argv[1])
        records = jsonReader.records()
        for count, record in enumerate(records, start=1):
            pass
        elapsed = time() - start
        print('Read %d JSON BLAST records in %.3f secs (%.0f records/sec)' % (
            count, elapsed, float(count) / float(elapsed)))

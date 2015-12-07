#!/usr/bin/env python

"""
Read simplified XML BLAST records and report the elapsed time.
"""

from __future__ import print_function

from Bio.Blast import NCBIXML
from time import time
import sys


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s file.json' % sys.argv[0], file=sys.stderr)
        sys.exit(1)
    else:
        start = time()
        with open(sys.argv[1]) as fp:
            for count, record in enumerate(NCBIXML.parse(fp)):
                pass
        elapsed = time() - start
        print('Read %d XML BLAST records in %.3f secs (%.0f records/sec)' % (
            count, elapsed, float(count) / float(elapsed)))

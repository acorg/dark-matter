#!/usr/bin/env python

from __future__ import print_function

from dark.blast import BlastRecords
from dark.blast.records import printBlastRecord


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: %s BLAST-hitfile...' % sys.argv[0], file=sys.stderr)
    else:
        blastRecords = BlastRecords(sys.argv[1:])
        for i, result in enumerate(blastRecords.records(), start=1):
            print('Read %d --------------------------' % i)
            printBlastRecord(result)

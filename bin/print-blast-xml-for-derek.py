#!/usr/bin/env python

from __future__ import print_function

from dark.blast import BlastRecords
from dark.blast.records import printBlastRecordForDerek

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: %s hitfiles...' % sys.argv[0], file=sys.stderr)
    else:
        blastRecords = BlastRecords(sys.argv[1:])
        hits = 0
        for i, result in enumerate(blastRecords.records(), start=1):
            oldHits = hits
            hits += printBlastRecordForDerek(result)

#!/usr/bin/env python

from dark.hacks import readBlastRecords, printBlastRecordForDerek

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Usage: %s hitfile' % sys.argv[0]
    else:
        results = readBlastRecords(sys.argv[1])
        hits = 0
        for i, result in enumerate(results, start=1):
            oldHits = hits
            hits += printBlastRecordForDerek(result)

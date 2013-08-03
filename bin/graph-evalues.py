#!/usr/bin/env python

from dark.utils import readBlastRecords, evalueGraph

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Usage: %s hitfile' % sys.argv[0]
    else:
        evalueGraph(readBlastRecords(sys.argv[1]), 10, 10,
                    find=lambda title: title.find('HKU4') > -1 and
                    title.find('complete genome') > -1, titles=False)

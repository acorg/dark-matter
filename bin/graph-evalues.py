#!/usr/bin/env python

"""
This script needs work!  Or will we never use it again?

There should be a --find option, or something.
"""

from __future__ import print_function

from dark.blast import BlastRecords
from dark.graphics import evalueGraph

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: %s BLAST-hitfile...' % sys.argv[0], file=sys.stderr)
    else:
        blastRecords = BlastRecords(sys.argv[1:])
        evalueGraph(blastRecords.records(), 10, 10,
                    find=lambda title: title.find('HKU4') > -1 and
                    title.find('complete genome') > -1, titles=False)

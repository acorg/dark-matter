#!/usr/bin/env python

from __future__ import print_function

import sys

from dark.blast import BlastRecords

if len(sys.argv) < 3:
    print('Usage: %s title record-file.(xml|json)...' % sys.argv[0],
          file=sys.stderr)
    sys.exit(1)
else:
    title = sys.argv[1]
    blastRecords = BlastRecords(sys.argv[2:])
    hits = blastRecords.filterHits(whitelist=set([title]),
                                   negativeTitleRegex='.')
    if title in hits.titles:
        print(hits.titles[title])
    else:
        print('Title %r was not hit by any reads.' % title)

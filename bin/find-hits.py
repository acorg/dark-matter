#!/usr/bin/env python

import sys

from dark.blast import BlastRecords

if len(sys.argv) < 3:
    print >>sys.stderr, (
        'Usage: %s title record-file.(xml|json)...' % sys.argv[0])
    sys.exit(1)
else:
    title = sys.argv[1]
    blastRecords = BlastRecords(sys.argv[2:])
    hits = blastRecords.filterHits(whitelist=set([title]),
                                   negativeTitleRegex='.')
    if title in hits.titles:
        print hits.titles[title]
    else:
        print 'Title %r was not hit by any reads.' % title

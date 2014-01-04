#!/usr/bin/env python

import sys

from dark.blast import BlastRecords

if len(sys.argv) != 3:
    print >>sys.stderr, (
        'Usage: %s record-file.(xml|json) title' % sys.argv[0])
    sys.exit(1)
else:
    recordFile, title = sys.argv[1:]
    blastRecords = BlastRecords(recordFile)
    hits = blastRecords.hits(whitelist=set([title]), negativeTitleRegex='.')
    if title in hits.titles:
        print hits.titles[title]
    else:
        print 'Title %r was not hit by any reads.' % title

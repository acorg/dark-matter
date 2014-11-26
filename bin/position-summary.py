#!/usr/bin/env python

from dark.summarize import summarizePosition
import sys


if len(sys.argv) != 3:
    print >>sys.stderr, 'Usage: %s file.fasta, position' % sys.argv[0]
    sys.exit(1)
else:
    fileName = sys.argv[1]
    position = int(sys.argv[2])

    result = summarizePosition(fileName, position)

    print 'There are %d sequences in %s' % (result['sequenceCount'], fileName)

    for base, count in result['countAtPosition'].items():
        print '%s: Total: %s; percentage: %s' % (
            base, count, (float(count) / float(result['sequenceCount'])))

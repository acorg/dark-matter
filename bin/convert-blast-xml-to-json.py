#!/usr/bin/env python

import sys
from os.path import basename

from dark.blast.conversion import XMLRecordsReader

if __name__ == '__main__':
    nArgs = len(sys.argv)
    if nArgs == 1:
        reader = XMLRecordsReader(sys.stdin)
        reader.saveAsJSON(sys.stdout)
    elif nArgs == 2:
        reader = XMLRecordsReader(sys.argv[1])
        reader.saveAsJSON(sys.stdout)
    elif nArgs == 3:
        reader = XMLRecordsReader(sys.argv[1])
        with open(sys.argv[2], 'w') as fp:
            reader.saveAsJSON(fp)
    else:
        print >>sys.stderr, (
            'Usage: %s [infile.xml [outfile.json]]' % basename(sys.argv[0]))
        sys.exit(1)

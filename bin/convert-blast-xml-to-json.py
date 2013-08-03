#!/usr/bin/env python

from dark.conversion import convertBlastXMLToJSON


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print >>sys.stderr, 'Usage: %s infile.xml' % sys.argv[0]
        sys.exit(1)
    else:
        convertBlastXMLToJSON(sys.argv[1], sys.stdout)

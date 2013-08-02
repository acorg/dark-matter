#!/usr/bin/env python

from dark.conversion import convertBlastXMLToJSON


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        blastFilename = sys.argv[1]
        jsonFilename = sys.stdout
    elif len(sys.argv) == 3:
        blastFilename = sys.argv[1]
        jsonFilename = sys.argv[2]
    else:
        print >>sys.stderr, 'Usage: %s infile.xml [outfile.json]' % sys.argv[0]
        sys.exit(1)

    convertBlastXMLToJSON(blastFilename, jsonFilename)

#!/usr/bin/env python

from dark.conversion import convertBlastXMLToJSON


def _convert(xmlfile, outfp):
    convertBlastXMLToJSON(xmlfile, outfp)


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        _convert(sys.argv[1], sys.stdout)
    if len(sys.argv) == 3:
        with open(sys.argv[2], 'w') as fp:
            _convert(sys.argv[1], fp)
    else:
        print >>sys.stderr, 'Usage: %s infile.xml [outfile.json]' % sys.argv[0]
        sys.exit(1)

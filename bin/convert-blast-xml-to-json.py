#!/usr/bin/env python

if __name__ == '__main__':
    from dark.conversion import XMLRecordsReader
    import sys
    if len(sys.argv) == 2:
        reader = XMLRecordsReader(sys.argv[1])
        reader.saveAsJSON(sys.stdout)
    if len(sys.argv) == 3:
        reader = XMLRecordsReader(sys.argv[1])
        with open(sys.argv[2], 'w') as fp:
            reader.saveAsJSON(sys.stdout)
    else:
        print >>sys.stderr, 'Usage: %s infile.xml [outfile.json]' % sys.argv[0]
        sys.exit(1)

#!/usr/bin/env python

import argparse
import bz2file
import sys

from dark.blast.conversion import XMLRecordsReader


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert a BLAST XML file to JSON.',
        epilog='Give a BLAST XML file and convert it to JSON. Optionally '
        'compress the JOSN output.'
    )

    parser.add_argument(
        'json', metavar='BLAST-JSON-file', nargs='?',
        help='the JSON filename where contents of "xml" should be written to.')

    parser.add_argument(
        '--xml', metavar='BLAST-XML-file', default=sys.stdin,
        help='the XML file of BLAST output.')

    parser.add_argument(
        '--bzip2', default=False, action='store_true',
        help='If True, compress the json output using bzip2.')

    args = parser.parse_args()

    if args.bzip2:
        fp = bz2file.BZ2File(args.json or sys.stdout, 'w')
    else:
        fp = open(args.json, 'w') if args.json else sys.stdout

    reader = XMLRecordsReader(args.xml)
    reader.saveAsJSON(fp)
    fp.close()

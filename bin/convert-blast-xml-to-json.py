#!/usr/bin/env python

import argparse
import bz2
import sys

from dark.blast.conversion import XMLRecordsReader


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert a BLAST XML file to JSON.',
        epilog='Give a BLAST XML file and convert it to JSON. Optionally '
        'compress the JOSN output file.'
    )

    parser.add_argument(
        'xml', metavar='BLAST-XML-file', type=str, help='the XML file of '
        'BLAST output.')

    parser.add_argument(
        'json', metavar='BLAST-JSON-file', type=str, nargs='?',
        help='the JSON filename where contents of "xml" should be written to.')

    parser.add_argument(
        '--compress', default=False, action='store_true',
        help='If True, compress the json output.')

    args = parser.parse_args()

    if args.compress:
        if args.json:
            reader = XMLRecordsReader(args.xml)
            fp = bz2.BZ2File(args.json, 'w')
            reader.saveAsJSON(fp)
            fp.close()
        else:
            raise ValueError('JSON filename must be specified if output '
                             'should be compressed.')
    else:
        if args.json:
            reader = XMLRecordsReader(args.xml)
            with open(args.json, 'w') as fp:
                reader.saveAsJSON(fp)
        else:
            reader = XMLRecordsReader(args.xml)
            reader.saveAsJSON(sys.stdout)

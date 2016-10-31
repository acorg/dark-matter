#!/usr/bin/env python

import argparse
import bz2file
import sys

from dark.diamond.conversion import DiamondTabularFormatReader


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert DIAMOND tabular output to JSON.',
        epilog=('Give a DIAMOND tabular file and convert it to JSON, '
                'optionally compressing the output. You *must* invoke '
                'DIAMOND with the following output specification: '
                '--outfmt 6 qtitle stitle bitscore evalue qframe qseq '
                'qstart qend sseq sstart send slen btop'
                'Note that only each line of the output is JSON (the full '
                'output is not valid JSON by itself).')
    )

    parser.add_argument(
        '--json', metavar='JSON-output-file',
        help=('The JSON filename to write the converted DIAMOND file to. If '
              'omitted, standard output will be used.'))

    parser.add_argument(
        '--diamond', metavar='DIAMOND-file', default=sys.stdin,
        help=('The DIAMOND tabular output file to convert. If omitted, '
              'standard input will be read.'))

    parser.add_argument(
        '--bzip2', default=False, action='store_true',
        help='Compress output using bzip2.')

    args = parser.parse_args()

    if args.bzip2:
        fp = bz2file.BZ2File(args.json or sys.stdout, 'w')
    else:
        fp = open(args.json, 'w') if args.json else sys.stdout

    reader = DiamondTabularFormatReader(args.diamond)
    reader.saveAsJSON(fp, writeBytes=args.bzip2)
    fp.close()

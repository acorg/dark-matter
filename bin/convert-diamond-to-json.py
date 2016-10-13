#!/usr/bin/env python

import argparse
import bz2file
import sys

from dark.blast.conversion import DiamondTabularFormatReader


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert a DIAMOND tabular output file to JSON.',
        epilog=('Give a DIAMOND tabular file and convert it to JSON. '
                'Optionally compress the JSON output. Make sure you have the '
                'right output flags set: The flags are: qseqid, sseqid, '
                'bitscore, evalue, qframe, qseq, qstart, qend, sseq, sstart, '
                'send, slen.')
    )

    parser.add_argument(
        '--json', metavar='JSON-output-file',
        help=('the JSON filename to write the converted DIAMOND file to. If '
              'omitted, standard output will be used.'))

    parser.add_argument(
        '--diamond', metavar='DIAMOND-file', default=sys.stdin,
        help=('the DIAMOND tabular output file to convert. If omitted, '
              'standard input will be read.'))

    parser.add_argument(
        '--bzip2', default=False, action='store_true',
        help='If True, compress the JSON output using bzip2.')

    args = parser.parse_args()

    if args.bzip2:
        fp = bz2file.BZ2File(args.json or sys.stdout, 'w')
    else:
        fp = open(args.json, 'w') if args.json else sys.stdout

    reader = DiamondTabularFormatReader(args.diamond)
    reader.saveAsJSON(fp)
    fp.close()

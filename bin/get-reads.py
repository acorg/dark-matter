#!/usr/bin/env python

import sys
from re import compile
import argparse

from dark.blast import BlastRecords


def main(recordFilenames, fastaFilename, title, xRange, eRange):
    """
    Prints the reads that are at a specified offset with a specified evalue.
    recordFilename: the result of a blast run, using outfmt 5.
    fastaFilename: the fastafile that was originally blasted.
    title: the title of the subject sequence, as output by BLAST.
    ranges: The first parameter must be a number of an interval on the
        x-axis from where the reads should be searched. The second parameter
        is optional and should be a converted value or an interval of
        converted evalues.
    """
    blastRecords = BlastRecords(recordFilenames, fastaFilename)
    hits = blastRecords.filterHits(whitelist=set([title]),
                                   negativeTitleRegex='.')
    if title not in hits.titles:
        print '%s: Title %r not found in BLAST output' % (sys.argv[0], title)
        sys.exit(3)

    hits.computePlotInfo()

    items = hits.titles[title]['plotInfo']['items']

    for item in items:
        hsp = item['hsp']
        if ((xRange is None or (xRange[0][0] <= hsp['subjectEnd'] and
                                xRange[0][1] >= hsp['subjectStart'])) and
                (eRange is None or
                    (eRange[0][0] <= item['convertedE'] <= eRange[0][1]))):
            print ('query: ', hits.fasta[item['readNum']].id, 'start: ',
                   hsp['subjectStart'], 'end: ', hsp['subjectEnd'],
                   'E-value: ', item['convertedE'])

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print >>sys.stderr, (
            'Usage: %s recordFilename, fastaFilename, '
            'title, xCoords, eCoords' % sys.argv[0])
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description=('Print the reads that are '
                     'at specified positions in an alignmentGraph'),
        epilog=('Given a JSON BLAST output file, a title and an x and / or '
                'eRange, print the reads that are within the given Ranges.'))

    parser.add_argument(
        'json', metavar='BLAST-JSON-file', type=str, nargs='+',
        help='the JSON file of BLAST output.')

    parser.add_argument(
        'fasta', metavar='fastaFile', type=str,
        help='the FASTA file of BLAST input.')

    parser.add_argument(
        'title', metavar='SEQUENCE-TITLE', type=str,
        help='The title of the subject sequence.')

    parser.add_argument(
        '--xRange', default=None,
        help='a range on the x-axis.')

    parser.add_argument(
        '--eRange', default=None,
        help='a range on the y-axis.')

    args = parser.parse_args()

    def _getRange(inputRange):
        if inputRange is None:
            return None
        rangeRegex = compile(r'^(\d+)(?:-(\d+))?$')
        ranges = []
        match = rangeRegex.match(inputRange)
        if match:
            start, end = match.groups()
            start = int(start)
            if end is None:
                end = start
            else:
                end = int(end)
            if start > end:
                start, end = end, start
            ranges.append((start, end))
        else:
            print >>sys.stderr, (
                'Illegal argument %r. Ranges must single numbers or '
                'number-number.' % inputRange)
            sys.exit(2)
        return ranges

    main(args.json, args.fasta, args.title,
         _getRange(args.xRange), _getRange(args.eRange))

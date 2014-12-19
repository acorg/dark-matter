#!/usr/bin/env python

import sys
from re import compile
import argparse

from dark.blast.alignments import BlastReadsAlignments
from dark.titles import TitlesAlignments
from dark.fasta import FastaReads


def main(recordFilename, fastaFilename, title, xRange, bitRange):
    """
    Prints the reads that are at a specified offset with a specified evalue.
    recordFilename: the result of a blast run, using outfmt 5.
    fastaFilename: the fastafile that was originally blasted.
    title: the title of the subject sequence, as output by BLAST.
    ranges: The first parameter must be a number of an interval on the
        x-axis from where the reads should be searched. The second parameter
        is optional and should be a converted value or an interval of
        converted bit scores.
    """
    reads = FastaReads(fastaFilename)
    blastReadsAlignments = BlastReadsAlignments(reads, recordFilename)
    filtered = blastReadsAlignments.filter(whitelist=set([title]),
                                           negativeTitleRegex='.')
    titlesAlignments = TitlesAlignments(filtered)

    if title not in titlesAlignments:
        print '%s: Title %r not found in BLAST output' % (sys.argv[0], title)
        sys.exit(3)

    for titleAlignment in titlesAlignments[title]:
        for hsp in titleAlignment.hsps:
            if ((xRange is None or (xRange[0][0] <= hsp.subjectEnd and
                                    xRange[0][1] >= hsp.subjectStart)) and
                (bitRange is None or (bitRange[0][0] <= hsp.score.score <=
                                      bitRange[0][1]))):
                print ('query: %s, start: %d, end: %d, score: %d' % (
                       titleAlignment.read.id, hsp.subjectStart,
                       hsp.subjectEnd, hsp.score.score))

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print >>sys.stderr, (
            'Usage: %s recordFilename, fastaFilename, '
            'title, xCoords, bitCoords' % sys.argv[0])
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description=('Print the reads that are '
                     'at specified positions in an alignmentGraph'),
        epilog=('Given a JSON BLAST output file, a title and an x and / or '
                'bitRange, print the reads that are within the given Ranges.'))

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
        '--bitRange', default=None,
        help='a bit score range on the y-axis.')

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
         _getRange(args.xRange), _getRange(args.bitRange))

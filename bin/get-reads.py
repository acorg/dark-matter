#!/usr/bin/env python

from __future__ import print_function

import sys
from re import compile
import argparse

from dark.blast.alignments import BlastReadsAlignments
from dark.titles import TitlesAlignments
from dark.fasta import FastaReads


def main(recordFilenames, fastaFilename, title, xRange, bitRange):
    """
    Print reads that match in a specified X-axis and bit score range.

    @param recordFilenames: A C{list} of C{str} file names contain results of a
        BLAST run, in JSON format.
    @param fastaFilename: The C{str} name of the FASTA file that was originally
        BLASTed.
    @param title: The C{str} title of the subject sequence, as output by BLAST.
    @param xRange: A (start, end) list of C{int}s, giving an X-axis range or
        C{None} if the entire X axis range should be printed.
    @param bitRange: A (start, end) list of C{int}s, giving a bit score range
        or C{None} if the entire bit score range should be printed.
    """
    reads = FastaReads(fastaFilename)
    blastReadsAlignments = BlastReadsAlignments(reads, recordFilenames)
    filtered = blastReadsAlignments.filter(whitelist=set([title]),
                                           negativeTitleRegex='.')
    titlesAlignments = TitlesAlignments(filtered)

    if title not in titlesAlignments:
        print('%s: Title %r not found in BLAST output' % (sys.argv[0], title))
        sys.exit(3)

    for titleAlignment in titlesAlignments[title]:
        for hsp in titleAlignment.hsps:
            if ((xRange is None or (xRange[0] <= hsp.subjectEnd and
                                    xRange[1] >= hsp.subjectStart)) and
                (bitRange is None or (bitRange[0] <= hsp.score.score <=
                                      bitRange[1]))):
                print(('query: %s, start: %d, end: %d, score: %d' % (
                       titleAlignment.read.id, hsp.subjectStart,
                       hsp.subjectEnd, hsp.score.score)))

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print((
            'Usage: %s recordFilename, fastaFilename, '
            'title, xCoords, bitCoords' % sys.argv[0]), file=sys.stderr)
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description=('Print the reads that are '
                     'at specified positions in an alignmentGraph'),
        epilog=('Given a JSON BLAST output file, a title and an x and / or '
                'bitRange, print the reads that are within the given Ranges.'))

    parser.add_argument(
        'json', metavar='BLAST-JSON-file', nargs='+',
        help='the JSON file of BLAST output.')

    parser.add_argument(
        'fasta', metavar='fastaFile', help='the FASTA file of BLAST input.')

    parser.add_argument(
        'title', metavar='SEQUENCE-TITLE',
        help='The title of the subject sequence.')

    parser.add_argument(
        '--xRange', default=None,
        help='a range on the x-axis.')

    parser.add_argument(
        '--bitRange', default=None,
        help='a bit score range on the y-axis.')

    args = parser.parse_args()

    def _getRange(inputRange):
        """
        Convert a string input range into a pair of integers.

        @param inputRange: A C{str}, either  a single number like '10', or a
            hyphen-separated pair of numbers like '3-9'. May also be C{None}.
        @return: Either C{None} if inputRange is C{None} or the empty string,
            else a (start, end) list of C{int}s.
        """
        if inputRange:
            rangeRegex = compile(r'^(\d+)(?:-(\d+))?$')
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
                return start, end
            else:
                print((
                    'Illegal argument %r. Ranges must single numbers or '
                    'number-number.' % inputRange), file=sys.stderr)
                sys.exit(2)

    main(args.json, args.fasta, args.title,
         _getRange(args.xRange), _getRange(args.bitRange))

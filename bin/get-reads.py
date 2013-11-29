#!/usr/bin/env python

from dark import utils
import sys
from re import compile
import argparse


def main(recordFilename, fastaFilename, hitId, xRange, eRange):
    """
    Prints the reads that are at a specified offset with a specified evalue.
    recordFilename: the result of a blast run, using outfmt 5.
    fastaFilename: the fasta file that was blasted.
    hitId: the hitId from the subject where the reads should be searched for.
        Must be of the format: 'gi|1234|abc|1234'
    ranges: The first parameter must be a number of an interval on the
        x-axis from where the reads should be searched. The second parameter
        is optional and should be a converted value or an interval of
        converted evalues.
    """
    allhits = utils.findHits(recordFilename, set([hitId]))

    fasta, summary = utils.summarizeHits(allhits, fastaFilename)

    hitInfo = summary[hitId]['items']
    if len(hitInfo) == 0:
        print >> sys.stderr, "%s: files are empty" % sys.argv[0]
        sys.exit(3)
    else:
        for item in hitInfo:
            hsp = item['hsp']
            if ((xRange is None or (xRange[0][0] <= hsp['subjectEnd'] and
                                    xRange[0][1] >= hsp['subjectStart'])) and
                    (eRange is None or
                        (eRange[0][0] <= item['convertedE'] <= eRange[0][1]))):
                print 'query: ', item['query'], 'start: ', hsp['subjectStart'],
                'end: ', hsp['subjectEnd'], 'E-value: ', item['convertedE']


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print >>sys.stderr, ('Usage: %s recordFilename, fastaFilename, \
            gi-number, xCoords, eCoords' % sys.argv[0])
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description='Print the reads that are \
        at specified positions in an alignmentGraph',
        epilog='Given a JSON BLAST output file, a FASTA sequence file, a hitId'
        'and a an x and / or eRange, print the reads that \
        are within the given Ranges.'
    )

    parser.add_argument(
        'json', metavar='BLAST-JSON-file', type=str,
        help='the JSON file of BLAST output.')

    parser.add_argument(
        'fasta', metavar='FASTA-file', type=str,
        help='the FASTA file of sequences that were given to BLAST.')

    parser.add_argument(
        'hitId', metavar='gi|1234|db|1234', type=str,
        help='The identifier of the subject.')

    parser.add_argument(
        '--xRange', default=None,
        help='a range on the x-axis.')

    parser.add_argument(
        '--eRange', default=None,
        help='a range on the y-axis.')

    args = parser.parse_args()

    def _getRange(inputRange):
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
                'number-number.' % arg)
            sys.exit(2)
        return ranges

    if args.xRange:
        xRange = _getRange(args.xRange)
    else:
        xRange = None

    if args.eRange:
        eRange = _getRange(args.eRange)
    else:
        eRange = None

    main(args.json, args.fasta, args.hitId, xRange, eRange)

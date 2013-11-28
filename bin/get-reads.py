#!/usr/bin/env python

from dark import utils
import sys
from re import compile


def main(recordFilename, fastaFilename, hitId, ranges):
    """
    Prints the reads that are at a specified offset with a specified evalue.
    recordFilename: the result of a blastrun, using outfmt 5.
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

    if hitInfo is None:
        print "files are empty"
        sys.exit(3)
    else:
        result = []
        evalueResult = []
        if ranges[0]:
            for item in hitInfo:
                if (ranges[0][0] <= item['hsp']['subjectEnd'] and
                        ranges[0][1] >= item['hsp']['subjectStart']):
                    result.append((item['query'], item['hsp']['subjectStart'],
                                  item['hsp']['subjectEnd'],
                                  item['convertedE']))

        try:
            if ranges[1]:
                for item in result:
                    if item[3] > ranges[1][0] and item[3] < ranges[1][1]:
                        evalueResult.append(item)
        except IndexError:
            pass

    if len(evalueResult) != 0:
        for item in evalueResult:
            print item
    else:
        for item in result:
            print item


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print >>sys.stderr, ('Usage: %s gi-number [offset1, offset2, ...]' %
                             sys.argv[0])
        sys.exit(1)

    rangeRegex = compile(r'^(\d+)(?:-(\d+))?$')
    ranges = []
    for arg in sys.argv[4:]:
        match = rangeRegex.match(arg)
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

    main(sys.argv[1], sys.argv[2], sys.argv[3], ranges)

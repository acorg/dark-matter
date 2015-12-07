#!/usr/bin/env python

from __future__ import print_function

from re import compile
import sys
from dark.entrez import getSequence


def main(gi, ranges):
    """
    Print the features of the genbank entry given by gi. If ranges is
    non-emtpy, only print features that include the ranges.

    gi: either a hit from a BLAST record, in the form
        'gi|63148399|gb|DQ011818.1|' or a gi number (63148399 in this example).
    ranges: a possibly empty list of ranges to print information for. Each
        range is a non-descending (start, end) pair of integers.
    """
    # TODO: Make it so we can pass a 'db' argument to getSequence.
    record = getSequence(gi)

    if record is None:
        print("Looks like you're offline.")
        sys.exit(3)
    else:
        printed = set()
        if ranges:
            for (start, end) in ranges:
                for index, feature in enumerate(record.features):
                    if (start < int(feature.location.end) and
                            end > int(feature.location.start) and
                            index not in printed):
                        print(feature)
                        printed.add(index)
        else:
            # Print all features.
            for feature in record.features:
                print(feature)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s gi-number [offset1, offset2, ...]' %
              sys.argv[0], file=sys.stderr)
        sys.exit(1)

    rangeRegex = compile(r'^(\d+)(?:-(\d+))?$')
    ranges = []
    for arg in sys.argv[2:]:
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
            print((
                'Illegal argument %r. Ranges must single numbers or '
                'number-number.' % arg), file=sys.stderr)
            sys.exit(2)

    main(sys.argv[1], ranges)

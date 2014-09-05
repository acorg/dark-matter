#!/usr/bin/env python

"""
Print a (Levenshtein) distance matrix for a set of known adaptors
given on the command line.

For an example of usage,
see https://notebooks.antigenic-cartography.org/terry/emc-adaptors.html
"""

from dark.distance import levenshtein

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Print a (Levenshtein) distance matrix for a set of '
                     'known adaptors'))

    parser.add_argument(
        'adaptors', type=str, nargs='+', metavar='adaptor',
        help='the set of adaptors that were used in sequencing')

    args = parser.parse_args()

    adaptors = args.adaptors
    nAdaptors = len(adaptors)
    length = len(adaptors[0])
    spaces = ' ' * length

    for i in xrange(length):
        print spaces,
        for adaptor in adaptors:
            print adaptor[i],
        print

    for i in xrange(nAdaptors):
        print adaptors[i],
        for j in xrange(nAdaptors):
            if j < i:
                print ' ',
            else:
                print levenshtein(adaptors[i], adaptors[j]),
        print

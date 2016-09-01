#!/usr/bin/env python

"""
Print a (Levenshtein) distance matrix for a set of known adaptors
given on the command line.

For an example of usage,
see https://notebooks.antigenic-cartography.org/terry/emc-adaptors.html
"""

from __future__ import print_function

from dark.distance import levenshtein

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Print a (Levenshtein) distance matrix for a set of '
                     'known adaptors'))

    parser.add_argument(
        'adaptors', nargs='+', metavar='adaptor',
        help='the set of adaptors that were used in sequencing')

    args = parser.parse_args()

    adaptors = args.adaptors
    nAdaptors = len(adaptors)
    length = len(adaptors[0])
    spaces = ' ' * length

    for i in range(length):
        print(spaces, end=' ')
        for adaptor in adaptors:
            print(adaptor[i], end=' ')
        print()

    for i in range(nAdaptors):
        print(adaptors[i], end=' ')
        for j in range(nAdaptors):
            if j < i:
                print(' ', end=' ')
            else:
                print(levenshtein(adaptors[i], adaptors[j]), end=' ')
        print()

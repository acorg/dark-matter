#!/usr/bin/env python

import sys
from dark.utils import getSeqFromGenbank


def main(gi, offsets):
    """
    Print the features of the genbank entry given by gi. If offsets is
    non-emtpy, only print features that include the offsets.

    gi: either a hit from a BLAST record, in the form
        'gi|63148399|gb|DQ011818.1|' or a gi number (63148399 in this example).
    offsets: a list (possibly empty) of offsets to print information for.
    """
    record = getSeqFromGenbank(gi)

    if record is None:
        print "Looks like you're offline."
        sys.exit(3)
    else:
        printed = set()
        if offsets:
            for offset in offsets:
                for index, feature in enumerate(record.features):
                    if offset in feature.location and index not in printed:
                        print feature
                        printed.add(index)
        else:
            # Print all features.
            for feature in record.features:
                print feature


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print >>sys.stderr, ('Usage: %s gi-number [offset1, offset2, ...]' %
                             sys.argv[0])
        sys.exit(1)
    try:
        offsets = map(int, sys.argv[2:])
    except ValueError:
        print >>sys.stderr, 'Offsets must be numeric.'
        sys.exit(2)
    else:
        main(sys.argv[1], offsets)

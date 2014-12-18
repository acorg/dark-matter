#!/usr/bin/env python

import sys

from dark.aa import find
from dark.aa import CODONS

from dark.aa import AA_LETTERS, ALL_PROPERTIES, PROPERTY_NAMES


def findOrDie(s):
    """
    Look up an amino acid.

    @param s: A C{str} amino acid specifier. This may be a full name,
        a 3-letter abbreviation or a 1-letter abbreviation. Case is ignored.
    @return: An C{AminoAcid} instance, if one can be found. Else exit.
    """
    aa = find(s)
    if aa:
        return aa
    else:
        print >>sys.stderr, 'Unknown amino acid or codon: %s' % s
        print >>sys.stderr, 'Valid arguments are: %s.' % CODONS.keys()
        sys.exit(1)


if len(sys.argv) == 1:
    # Show info on all amino acides.
    aas = AA_LETTERS
else:
    aas = sys.argv[1:]

for aa in map(findOrDie, aas):
    print aa.name
    print '  3-letter abbreviation: %s' % aa.abbrev3
    print '  1-letter abbreviation: %s' % aa.abbrev1
    print '  Codons: %s' % ', '.join(sorted(aa.codons))

    properties = []
    print '  Properties:',
    for prop in ALL_PROPERTIES:
        if aa.properties & prop:
            properties.append(PROPERTY_NAMES[prop])
    print ', '.join(properties)

    print '  Property details:'
    for propertyDetail, value in aa.propertyDetails.items():
        print '    %s: %s' % (propertyDetail, value)

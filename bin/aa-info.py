#!/usr/bin/env python

import sys
from pprint import pprint

from dark.aa import find
from dark.aa import CODONS

from dark.aa import (
    ACIDIC, ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC,
    HYDROPHOBIC, HYDROXYLIC, NEGATIVE, POLAR, SMALL, SULPHUR, TINY)

INDIVIDUAL_PROPERTIES = {
    ACIDIC: 'Acidic',
    ALIPHATIC: 'Aliphatic',
    AROMATIC: 'Aromatic',
    BASIC_POSITIVE: 'Basic positive',
    HYDROPHILIC: 'Hydrophilic',
    HYDROPHOBIC: 'Hydrophobic',
    HYDROXYLIC: 'Hydroxylic',
    NEGATIVE: 'Negative',
    POLAR: 'Polar',
    SMALL: 'Small',
    SULPHUR: 'Sulphur',
    TINY: 'Tiny',
}


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
    print >>sys.stderr, 'Usage: %s amino_acid...' % sys.argv[0]
    sys.exit(1)

for aa in map(findOrDie, sys.argv[1:]):
    print 'Name: %s' % aa.name
    print 'Abbrev3: %s' % aa.abbrev3
    print 'Abbrev1: %s' % aa.abbrev1
    print 'Codons: %r' % aa.codons

    properties = []
    print 'Properties:',
    for property_, name in INDIVIDUAL_PROPERTIES.items():
        if aa.properties & property_:
            properties.append(name)
    print ', '.join(properties)

    print 'Property details:'
    pprint(aa.propertyDetails, indent=2)

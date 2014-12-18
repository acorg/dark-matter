#!/usr/bin/env python

import sys

from dark.aa import find
from dark.codonDistance import codonInformation
from dark.aa import CODONS


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


if len(sys.argv) != 3:
    print >>sys.stderr, 'Usage: %s amino_acid amino_acid' % sys.argv[0]
    sys.exit(1)
else:
    aa1 = findOrDie(sys.argv[1])
    aa2 = findOrDie(sys.argv[2])

    distances = codonInformation(aa1.codons, aa2.codons)
    print 'Codons:'
    print '%s: %s' % (aa1.name, ', '.join(sorted(aa1.codons)))
    print '%s: %s' % (aa2.name, ', '.join(sorted(aa2.codons)))

    print 'Distances:'
    sortedDistances = sorted(distances.keys())
    for distance in sortedDistances:
        for codon in distances[distance]:
            print '%d: %s: %s -> %s: %s' % (distance, aa1.name, codon[0],
                                            aa2.name, codon[1])

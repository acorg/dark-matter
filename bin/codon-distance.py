#!/usr/bin/env python

from dark.codonDistance import codonInformation
import sys


if len(sys.argv) != 3:
    print >>sys.stderr, 'Usage: %s amino_acid, amino_acid' % sys.argv[0]
    sys.exit(1)
else:
    aa1 = sys.argv[1]
    aa2 = sys.argv[2]

    result = codonInformation(aa1, aa2)

    print 'Codons:', aa1 + ':',
    for c in result['codons'][aa1]:
        print c,

    print aa2 + ':',
    for c in result['codons'][aa2]:
        print c,
    print '\n'
    print 'Minimum distance: %d, Codons: %s, %s' % (result['minDistance'],
                                                    result['minCodons'][0],
                                                    result['minCodons'][0])
    print 'Maximum distance: %d, Codons: %s, %s' % (result['maxDistance'],
                                                    result['maxCodons'][0],
                                                    result['maxCodons'][0])

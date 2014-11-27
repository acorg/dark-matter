#!/usr/bin/env python

import sys

from dark.codonDistance import codonInformation
from dark.aa import CODONS


if len(sys.argv) != 3:
    print >>sys.stderr, 'Usage: %s amino_acid amino_acid' % sys.argv[0]
    sys.exit(1)
else:
    aa1 = sys.argv[1]
    aa2 = sys.argv[2]

    try:
        codonsAA1 = CODONS[aa1.lower()]
        codonsAA2 = CODONS[aa2.lower()]
    except KeyError:
        print >>sys.stderr, ('Your arguments need to be two of: %s') % (
            CODONS.keys())
        sys.exit(1)

    result = codonInformation(codonsAA1, codonsAA2)
    print 'Codons:'
    print '%s: %s' % (aa1, ' '.join(codonsAA1))
    print '%s: %s' % (aa2, ' '.join(codonsAA2))

    print 'Distances:'
    for i in result:
        for codon in result[i]:
            print '%d: %s: %s; %s: %s' % (i, aa1, codon[0], aa2, codon[1])

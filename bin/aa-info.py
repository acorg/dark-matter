#!/usr/bin/env python

from __future__ import print_function

import sys

from dark.aa import find
from dark.aa import CODONS

from dark.aa import AA_LETTERS, ALL_PROPERTIES, PROPERTY_NAMES

args = AA_LETTERS if len(sys.argv) == 1 else sys.argv[1:]
error = False

for arg in args:
    aas = list(find(arg))
    if aas:
        for aa in aas:
            print(aa.name)
            print('  3-letter abbreviation: %s' % aa.abbrev3)
            print('  1-letter abbreviation: %s' % aa.abbrev1)
            print('  Codons: %s' % ', '.join(sorted(aa.codons)))

            properties = []
            print('  Properties:', end=' ')
            for prop in ALL_PROPERTIES:
                if aa.properties & prop:
                    properties.append(PROPERTY_NAMES[prop])
            print(', '.join(properties))

            print('  Property details:')
            for propertyDetail, value in aa.propertyDetails.items():
                print('    %s: %s' % (propertyDetail, value))
    else:
        error = True
        print('Unknown amino acid or codon: %s' % arg, file=sys.stderr)

if error:
    print('Valid arguments are: %s.' % list(CODONS), file=sys.stderr)

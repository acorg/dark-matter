#!/usr/bin/env python

import sys

from dark.aa import find
from dark.aaVars import CODONS
from dark.codonDistance import codonInformation


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
        print("Unknown amino acid or codon: %s" % s, file=sys.stderr)
        print("Valid arguments are: %s." % ", ".join(sorted(CODONS)), file=sys.stderr)
        sys.exit(1)


def codonKey(codonInfo):
    return codonInfo[0], codonInfo[1]


def transitionsKey(codonInfo):
    return codonInfo[2]


if len(sys.argv) != 3:
    print("Usage: %s amino_acid amino_acid" % sys.argv[0], file=sys.stderr)
    sys.exit(1)
else:
    (aa1,) = findOrDie(sys.argv[1])
    (aa2,) = findOrDie(sys.argv[2])

    distances = codonInformation(aa1.codons, aa2.codons, countTransitions=True)
    print("Codons:")
    print("%s: %s" % (aa1.name, ", ".join(sorted(aa1.codons))))
    print("%s: %s" % (aa2.name, ", ".join(sorted(aa2.codons))))

    # Output is sorted first by distance, then by number of transitions
    # (highest to lowest) then by codon DNA. The idea being to present the
    # possible codon changes to get from one amino acid to another in the
    # order that requires the least change to the most.

    print("Distances:")
    for distance in sorted(distances):
        for codon1, codon2, transitions in sorted(
            sorted(distances[distance], key=codonKey), key=transitionsKey, reverse=True
        ):
            print(
                f"{distance}: {aa1.name}: {codon1} -> "
                f"{aa2.name}: {codon2} (transitions: {transitions}, "
                f"transversions {distance - transitions})"
            )

#!/usr/bin/env python

import sys

from rich import print as rprint

from dark.aa import AminoAcid, find
from dark.aaVars import CODONS, PROPERTY_NAMES, STOP_CODONS
from dark.codonDistance import codonInformation

# Nucleotide colors. Note that 'A' will just be rendered as bold text in
# whatever the default terminal color is (often black or white, depending on
# the user's preferences). This is not perfect, but it's hard to be perfect.
COLORS = {
    "A": "",
    "C": "red",
    "G": "green",
    "T": "blue",
}


def colorCodons(codon1: str, codon2: str) -> tuple[str, str]:
    """
    Take two codons and return the same two codons but with differing
    nucleotides colored.
    """
    result1 = result2 = ""

    for nt1, nt2 in zip(codon1, codon2):
        if nt1 == nt2:
            result1 += nt1
            result2 += nt2
        else:
            color1 = f"bold {COLORS[nt1]}"
            color2 = f"bold {COLORS[nt2]}"
            result1 += f"[{color1}]{nt1}[/{color1}]"
            result2 += f"[{color2}]{nt2}[/{color2}]"

    return result1, result2


def findOrDie(s) -> AminoAcid:
    """
    Look up an amino acid.

    @param s: A C{str} amino acid specifier. This may be a full name,
        a 3-letter abbreviation or a 1-letter abbreviation. Case is ignored.
    @return: An C{AminoAcid} instance, if one can be found. Else exit.
    """
    aas = list(find(s))
    if len(aas) == 1:
        return aas[0]

    if aas:
        print(
            f"Multiple amino acids matched {s!r}:",
            ", ".join(sorted(aa.name for aa in aas)),
            file=sys.stderr,
        )
    else:
        if s.upper() in STOP_CODONS:
            print(f"Argument {s!r} is a stop codon!", file=sys.stderr)
        else:
            print(f"Unknown amino acid or codon: {s!r}", file=sys.stderr)
        print(
            "Valid amino acid arguments are: %s." % ", ".join(sorted(CODONS)),
            file=sys.stderr,
        )

    sys.exit(1)


def main():
    def codonKey(codonInfo):
        return codonInfo[0], codonInfo[1]

    def transitionsKey(codonInfo):
        return codonInfo[2]

    if len(sys.argv) != 3:
        print("Usage: %s amino_acid amino_acid" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    arg1 = sys.argv[1]
    arg2 = sys.argv[2]

    if arg1 == arg2:
        print("Identical arguments, nothing to do.", file=sys.stderr)
        sys.exit(1)

    aa1 = findOrDie(arg1)
    aa2 = findOrDie(arg2)

    if aa1.name == aa2.name:
        print(
            f"Arguments {arg1!r} and {arg2!r} both correspond to {aa1.name}, "
            "nothing to do.",
            file=sys.stderr,
        )
        sys.exit(1)

    print("Names:")
    print(f"  {arg1}: {aa1.name}")
    print(f"  {arg2}: {aa2.name}")

    def propertyNames(properties: set[int]) -> str:
        return ", ".join(PROPERTY_NAMES[p] for p in properties)

    aa1Properties = aa1.propertySet()
    aa2Properties = aa2.propertySet()

    print("Properties:")
    print(f"  {aa1.name}: {propertyNames(aa1Properties)}")
    print(f"  {aa2.name}: {propertyNames(aa2Properties)}")

    # aa1Only = propertyNames(aa1Properties - aa2Properties)
    # if aa1Only:
    #     print(f"  {aa1.name} only: {aa1Only}")

    # aa2Only = propertyNames(aa2Properties - aa1Properties)
    # if aa2Only:
    #     print(f"  {aa2.name} only: {aa2Only}")

    common = propertyNames(aa1Properties & aa2Properties) or "Nothing"
    print(f"  In common: {common}")

    print("Codons:")
    print("  %s: %s" % (aa1.name, ", ".join(sorted(aa1.codons))))
    print("  %s: %s" % (aa2.name, ", ".join(sorted(aa2.codons))))

    # Codon distance output is sorted first by distance, then by number of
    # transitions (lowest to highest) then by codon DNA. The idea being to
    # present the possible codon changes to get from one amino acid to
    # another in the order that requires the least change to the most.

    print("Codon distances:")
    distances = codonInformation(aa1.codons, aa2.codons, countTransitions=True)

    for distance in sorted(distances):
        for codon1, codon2, transitions in sorted(
            sorted(distances[distance], key=codonKey), key=transitionsKey, reverse=True
        ):
            transversions = distance - transitions
            codon1, codon2 = colorCodons(codon1, codon2)
            rprint(
                f"  {distance}: {aa1.name}: {codon1} -> {aa2.name}: {codon2}",
                end=" ",
            )
            print(f"(transitions: {transitions}, transversions {transversions})")


if __name__ == "__main__":
    main()

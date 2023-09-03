#!/usr/bin/env python

import sys
from argparse import ArgumentParser

from dark.aa import find
from dark.aaVars import (
    CODONS,
    NAMES_TO_ABBREV1,
    ALL_PROPERTIES,
    PROPERTY_NAMES,
    STOP_CODONS,
)

parser = ArgumentParser()
parser.add_argument("targets", nargs="*", help="The items to look up.")
parser.add_argument(
    "--details", action="store_true", help="Display numeric amino acid property values."
)
args = parser.parse_args()

targets = args.targets or sorted(NAMES_TO_ABBREV1)

error = False

for target in targets:
    aas = list(find(target))
    if aas:
        for aa in aas:
            print(aa.name)
            print("  3-letter abbreviation: %s" % aa.abbrev3)
            print("  1-letter abbreviation: %s" % aa.abbrev1)
            print("  Codons: %s" % ", ".join(sorted(aa.codons)))

            properties = []
            print("  Properties:", end=" ")
            for prop in ALL_PROPERTIES:
                if aa.properties & prop:
                    properties.append(PROPERTY_NAMES[prop])
            print(", ".join(properties))

            if args.details:
                print("  Property details:")
                for propertyDetail, value in aa.propertyDetails.items():
                    print("    %s: %s" % (propertyDetail, value))
    else:
        if (
            target.upper() in STOP_CODONS
            or target.upper().replace("U", "T") in STOP_CODONS
        ):
            print("Stop codon.")
        else:
            error = True
            print("Unknown amino acid or codon: %s" % target, file=sys.stderr)

if error:
    print("Valid AA arguments are: %s." % list(CODONS), file=sys.stderr)
    sys.exit(1)

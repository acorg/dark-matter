#!/usr/bin/env python

import argparse

from Bio.Seq import Seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Print the amino acid translation of one or more DNA sequences."
    )
    parser.add_argument(
        "sequences",
        metavar="DNA-SEQUENCE",
        nargs="+",
        help="A DNA sequence to translate.",
    )
    args = parser.parse_args()

    for dna in args.sequences:
        print(Seq(dna).translate())

#!/usr/bin/env python

import argparse
import sys

from dark.analyze_reads import getPrefixAndSuffix, trimReads
from dark.reads import Reads

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Given a FASTA file, detect a prefix/suffix common to all reads "
            "(e.g., leftover primer or adaptor sequence) and print the reads "
            "with that trimmed off."
        )
    )
    parser.add_argument("filename", help="The FASTA file to trim.")
    args = parser.parse_args()

    filename = args.filename
    with open(filename, "rb") as fp:
        prefix, suffix = getPrefixAndSuffix(fp)

    print(f"Prefix length {prefix}, suffix length {suffix}", file=sys.stderr)

    with open(filename, "rb") as fp:
        Reads(trimReads(prefix, suffix, fp)).save(sys.stdout)

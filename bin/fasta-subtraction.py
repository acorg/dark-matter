#!/usr/bin/env python

import argparse

from dark.fasta import fastaSubtract

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Given two or more FASTA file names as arguments, remove from "
            "the first file any read whose id also appears in a later "
            "file, and print what's left."
        )
    )

    parser.add_argument(
        "fastaFiles",
        nargs="+",
        metavar="FASTA-file",
        help=(
            "A FASTA file name. Reads are removed from the first file "
            "named if their id also appears in any of the later files."
        ),
    )

    args = parser.parse_args()

    for read in fastaSubtract(args.fastaFiles):
        print(read.toString(), end="")

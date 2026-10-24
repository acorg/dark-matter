#!/usr/bin/env python

import argparse

from dark.reads import DNARead

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Print the reverse complement of one or more DNA sequences."
    )
    parser.add_argument(
        "sequences", metavar="SEQUENCE", nargs="+", help="A DNA sequence."
    )
    args = parser.parse_args()

    for seq in args.sequences:
        print(DNARead("id", seq).reverseComplement().sequence)

#!/usr/bin/env python

"""
Given a set of FASTA sequence identifiers from sys.argv and/or in a file, read
FASTA from stdin, and print FASTA to stdout for the given sequence ids.
"""

import argparse
import sys

from dark.fasta import FastaReads
from dark.reads import Reads


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract a subset of FASTA reads by id",
        epilog=(
            "Given a set of FASTA sequence identifiers from sys.argv "
            "or in a file, read FASTA from stdin, and print FASTA to stdout "
            "for the given sequence ids."
        ),
    )

    parser.add_argument("ids", default=None, nargs="*", help="Wanted read ids.")

    parser.add_argument("--readIdFile", default=None, help="A file of wanted read ids.")

    args = parser.parse_args()
    wanted = set()

    if args.ids:
        wanted.update(args.ids)

    if args.readIdFile is not None:
        with open(args.readIdFile) as fp:
            for line in fp:
                wanted.add(line[:-1])

    found = []

    if wanted:
        for seq in FastaReads(sys.stdin):
            if seq.id in wanted:
                wanted.remove(seq.id)
                found.append(seq)

        if found:
            Reads(found).save(sys.stdout)

        print("Found %d sequences." % len(found), file=sys.stderr)

        if wanted:
            print(
                "WARNING: %d sequence%s not found: %s"
                % (len(wanted), "" if len(wanted) == 1 else "s were", " ".join(wanted)),
                file=sys.stderr,
            )
    else:
        # No wanted ids were given.
        parser.print_help(sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

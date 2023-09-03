#!/usr/bin/env python

from random import choice
from operator import itemgetter

from dark.aaVars import CODONS
from dark.reads import Read, addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given amino acid FASTA on stdin, write FASTA with DNA "
            "sequences to stdout, where the DNA sequences would "
            "translate to the given amino acid seqeunces. The DNA "
            "bases used for each amino acid is either fixed or is "
            "chosen at random from the corresponding codons."
        )
    )

    parser.add_argument(
        "--random",
        action="store_true",
        help=(
            "Use codons chosen at random for each amino acid. Otherwise a "
            "fixed codon will be used for each amino acid so that the "
            "output for a given amino acid sequence is always identical."
        ),
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    getCodon = choice if args.random else itemgetter(0)

    for aaRead in reads:
        dnaRead = Read(
            aaRead.id, "".join(getCodon(CODONS[aa]) for aa in aaRead.sequence)
        )
        print(dnaRead.toString("fasta"), end="")

#!/usr/bin/env python

from Bio.Seq import translate

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Translate FASTA and write the result to standard output."
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()

    for read in parseFASTACommandLineOptions(args):
        print(f">{read.id}\n{translate(read.sequence)}")

#!/usr/bin/env python

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin, write the sequence ids and lengths to stdout."
        )
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        print("%s %d" % (read.id, len(read)))

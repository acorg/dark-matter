#!/usr/bin/env python

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


import argparse

parser = argparse.ArgumentParser(
    description=(
        "Given FASTA/Q on stdin, write FASTQ to stdout with all nucleotide quality "
        "scores set to a given value."
    )
)

parser.add_argument(
    "--quality",
    "-q",
    type=int,
    default=30,
    help="The integer PHRED quality to set. Sanger (32-based PHRED is used).",
)

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)
quality = chr(ord(" ") + args.quality)

for read in reads:
    read.quality = quality * len(read)
    print(read.toString(format_="fastq"), end="")

#!/usr/bin/env python

import sys
import argparse

from dark.reads import (
    addFASTACommandLineOptions,
    parseFASTACommandLineOptions,
    simpleReadSplitter,
)


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Split FASTA/Q sequences longer than a given length. Print the beginning "
        "and end sequences of sufficiently long sequences separately."
    ),
)

parser.add_argument(
    "--length",
    type=int,
    required=True,
    help="Sequences longer than this will be split in two.",
)

parser.add_argument(
    "--leftPrefix",
    default="",
    help="The string to put before the id of the starting (left) subsequence.",
)

parser.add_argument(
    "--leftSuffix",
    default="",
    help="The string to put after the id of the starting (left) subsequence.",
)

parser.add_argument(
    "--rightPrefix",
    default="",
    help="The string to put before the id of the ending (right) subsequence.",
)

parser.add_argument(
    "--rightSuffix",
    default="",
    help="The string to put after the id of the ending (right) subsequence.",
)

parser.add_argument(
    "--noWarn",
    dest="warn",
    action="store_false",
    help=(
        "Do not warn the user if no prefixes are given (which can result in "
        "FASTA output with duplicated IDs."
    ),
)

addFASTACommandLineOptions(parser)
args = parser.parse_args()

assert args.length > 0, "The --length option must be positive."

if args.warn and (
    (args.leftPrefix == "" and args.leftSuffix == "")
    or (args.rightPrefix == "" and args.rightSuffix == "")
):
    print(
        "WARNING: You have not specified enough left/right ID sequence fragment "
        "prefixes/suffixes. This means found fragments (if any) will receive ids "
        "that are identical to the longer sequence they are derived from.",
        file=sys.stderr,
    )

reads = parseFASTACommandLineOptions(args)

splitter = simpleReadSplitter(
    args.length, args.leftPrefix, args.leftSuffix, args.rightPrefix, args.rightSuffix
)

for read in parseFASTACommandLineOptions(args):
    for fragment in splitter(read):
        print(fragment.toString(), end="")

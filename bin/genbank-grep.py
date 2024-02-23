#!/usr/bin/env python

import sys
import re
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Read a GenBank flat file and print GenBank or FASTA for matching records."
    ),
)

parser.add_argument(
    "genbankFile", metavar="FILE.gb", nargs="+", help="The GenBank file(s) to examine."
)

parser.add_argument(
    "--idRegex",
    metavar="REGEX",
    default=".",
    help=(
        "Only print records whose id (or description, if --matchDescription "
        "is used) matches this regular expression."
    ),
)

parser.add_argument(
    "--match",
    choices=("fullmatch", "match", "search"),
    default="search",
    help="How to match the regular expression. See pydoc re for an " "explanation.",
)

parser.add_argument(
    "--matchDescription",
    action="store_true",
    help="Match on the record description as well as the id.",
)

parser.add_argument(
    "--invert",
    "-v",
    action="store_true",
    help="Invert the match. Only print records that do not match.",
)

parser.add_argument(
    "--ignoreCase", "-i", action="store_true", help="Use case-insensitive matching."
)

parser.add_argument(
    "--format", choices=("fasta", "gb"), default="gb", help="The output format."
)

parser.add_argument(
    "--expectedCount",
    type=int,
    help=(
        "The expected number of records to find. If this number is not "
        "found, print a message and exit non-zero."
    ),
)

parser.add_argument(
    "--verbose",
    action="store_true",
    help="Print a summary of records read and printed.",
)

args = parser.parse_args()

recordCount = matchCount = 0
match = getattr(
    re.compile(args.idRegex, flags=(re.I if args.ignoreCase else 0)), args.match
)

matchDescription = args.matchDescription
invert = args.invert


def filterFunc(records):
    global matchCount, recordCount
    for record in records:
        recordCount += 1
        id_ = record.id + (f" {record.description}" if record.description else "")
        if match(id_ if matchDescription else record.id):
            emit = not invert
        else:
            emit = invert

        if emit:
            matchCount += 1
            yield record


for filename in args.genbankFile:
    with open(filename) as fp:
        SeqIO.write(filterFunc(SeqIO.parse(fp, "gb")), sys.stdout, args.format)

if args.verbose:
    print(
        f'Read {recordCount} record{"s" if recordCount > 1 else ""}, '
        f"printed {matchCount}",
        file=sys.stderr,
    )

if args.expectedCount is not None and matchCount != args.expectedCount:
    print(
        f"Expected {args.expectedCount} records but found {matchCount}. " f"Exiting.",
        file=sys.stderr,
    )
    sys.exit(1)

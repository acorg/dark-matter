#!/usr/bin/env python

import sys
import argparse
import re

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions, Reads

parser = argparse.ArgumentParser(
    description=(
        "Write sorted FASTA/Q to stdout. Sorting is by sequence id, "
        "then by sequence, then by quality (if FASTQ)"
    )
)

# TODO: Support multiple regexes and apply them one by one to allow matching on
# different regions in an arbitrary order.
parser.add_argument(
    "--regex",
    help=(
        "A regular expression to specify a region (or regions) to extract "
        "from sequence ids to sort on. Regions will be converted to "
        "integers if possible."
    ),
)


parser.add_argument(
    "--reverse", "-r", action="store_true", help="Sort by decreasing value."
)

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)

if args.regex:
    regex = re.compile(args.regex)

    if regex.groups == 0:
        print(
            "You have passed a regular expression that has no capturing "
            "group(s) specified using (...).",
            file=sys.stderr,
        )
        sys.exit(1)

    def key(read):
        result = []
        match = regex.search(read.id)
        if match:
            for text in match.groups():
                try:
                    value = int(text)
                except ValueError:
                    value = text
                result.append(value)
        return result

else:

    def key(read):
        return read.id


Reads(sorted(reads, key=key, reverse=args.reverse)).save(
    sys.stdout, "fastq" if args.fastq else "fasta"
)

#!/usr/bin/env python

import sys
import argparse

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Find sequences that have a specific property (as detected "
        'by the "relabel" function in the file specified by '
        "--relabelFunctionFile) and print a TAB-separated list of "
        "renamings to standard output."
    ),
)

parser.add_argument(
    "--relabelFunctionFile",
    required=True,
    help=(
        'A file containing a Python function named "relabel" that takes a '
        "list of reads and returns a dictionary of renamings."
    ),
)

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)

relabel = origRelabel = object()

try:
    exec(open(args.relabelFunctionFile).read())
except Exception as e:
    print(
        "Could not execute Python in %s: %s" % (args.relabelFunctionFile, e),
        file=sys.stderr,
    )
    sys.exit(1)

if relabel is origRelabel:
    print(
        'The Python file %s did not define a function called "relabel".'
        % args.relabelFunctionFile,
        file=sys.stderr,
    )
    sys.exit(2)

for old, new in relabel(reads).items():
    print("%s\t%s" % (old, new))

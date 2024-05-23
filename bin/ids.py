#!/usr/bin/env python

import sys
import argparse

from dark.idutils import ids


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "spec",
        nargs="+",
        help=(
            "One or two arguments to specify the start and (optionally) the end "
            "range of ids to create."
        ),
        metavar="ID",
    )

    parser.add_argument(
        "--prefix",
        help=(
            "The id prefix. If this is not given and two arguments are supplied, "
            "their longest common prefix will be used."
        ),
    )

    parser.add_argument(
        "--sep", default=" ", help="The string to print between successive ids."
    )

    parser.add_argument(
        "-0",
        "--zeroes",
        action="store_true",
        help="Print numeric suffixes with leading zeroes.",
    )

    parser.add_argument(
        "-n",
        "--count",
        dest="maxResults",
        type=int,
        help="The number of items to generate.",
    )

    args = parser.parse_args()

    if len(args.spec) == 1:
        start, end = args.spec[0], None
    elif len(args.spec) == 2:
        start, end = args.spec
    else:
        print("I only know how to process one or two arguments.", file=sys.stderr)
        sys.exit(1)

    sep = {
        r"\t": "\t",
        r"\n": "\n",
    }.get(args.sep, args.sep)

    # Iterate through the yielded results instead of trying to compute them
    # all at once and joining the resulting list.
    first = True

    for id_ in ids(start, end, args.prefix, sep, args.zeroes, args.maxResults):
        if not first:
            print(sep, end="")
        print(f"{id_}", end="")
        first = False

    if not first:
        print()


if __name__ == "__main__":
    main()

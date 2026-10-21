#!/usr/bin/env python

import argparse
import sys

from dark.reads import (
    ReadsInRAM,
    addFASTACommandLineOptions,
    parseFASTACommandLineOptions,
)


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find nucleotide sequences bounding a region in an MSA and merge "
            "the gaps within that region. I.e., merge all the gaps into a single "
            "gap, place it at a specific offset within the region, and put the "
            "non-gap characters (now contiguous) to the left or right (or around) "
            "the gap."
        ),
    )

    addFASTACommandLineOptions(parser)

    parser.add_argument("--prefix", help="The start of the match to find in sequences.")

    parser.add_argument("--suffix", help="The end of the match to find in sequences.")

    parser.add_argument(
        "--id",
        help=(
            "The id of the input sequence in which the prefix and suffix must be "
            "found. If not specified, all input sequences will be examined and "
            "the offsets in which the prefix/suffix are found must all agree."
        ),
    )

    parser.add_argument(
        "--gapOffset",
        type=int,
        default=0,
        help=(
            "The offset at which the merged gaps should be put. The default (0) will "
            "put the merged gaps on the left (i.e., with zero non-gap chars to the "
            "left of them. This value cannot be larger than the number of non-gap "
            "chars in the sequence region (otherwise there are not enough non-gap "
            "chars to put before the gap). Use -1 to put the merged gaps on the right."
        ),
    )

    parser.add_argument(
        "--allowDistributedMatch",
        action="store_true",
        help=(
            "If --prefix and --suffix are both used, the default behaviour is to "
            "make sure that they are both found in at least one of the input "
            "sequences. If you don't mind if the prefix is found in one sequence "
            "and the suffix in another, use this option to allow it."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print a summary of prefixes/suffixes found to standard error.",
    )

    parser.add_argument(
        "--zeroBased",
        action="store_true",
        help=(
            "If --tsv or --verbose are used, the reported prefix and suffix positions "
            "will refer to zero-based Python style [closed, open) regions."
        ),
    )

    parser.add_argument(
        "--gapCharacter", default="-", help="The 1-letter alignment gap character."
    )

    parser.add_argument(
        "--tsv", help="The (optional) name of a TSV file to write match offsets to."
    )

    return parser.parse_args()


def main():
    args = get_args()
    reads = ReadsInRAM(parseFASTACommandLineOptions(args))
    result, summary = reads.mergeGapsInRegion(
        id_=args.id,
        prefix=args.prefix,
        suffix=args.suffix,
        tsv=args.tsv,
        allowDistributedMatch=args.allowDistributedMatch,
        gapCharacter=args.gapCharacter,
        gapOffset=args.gapOffset,
        adjust=int(not args.zeroBased),
        summarize=args.verbose,
    )

    result.save(sys.stdout)

    if args.verbose:
        assert summary is not None
        print(summary, file=sys.stderr)


if __name__ == "__main__":
    main()

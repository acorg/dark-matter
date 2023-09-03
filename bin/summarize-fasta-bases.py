#!/usr/bin/env python

from collections import defaultdict
from math import log10

from dark.aaVars import NAMES
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions
from dark.summarize import sequenceCategoryLengths


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin, write a summary of sequence base "
            "categories to stdout. It is currently not possible to "
            "specify the categories on the command line."
        )
    )

    parser.add_argument(
        "--baseType",
        default="nucl",
        choices=("nucl", "prot"),
        help="The type of the bases in the input.",
    )

    parser.add_argument(
        "--minLength",
        default=1,
        type=int,
        help=(
            "If specified, stretches of reads that are less than this "
            "length will not be reported but will be summarized by an "
            "ellipsis."
        ),
    )

    parser.add_argument(
        "--concise",
        action="store_true",
        default=False,
        help="If specified, do not show the individual sequence regions.",
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    if args.baseType == "nucl":
        categories = {
            "A": "nucl",
            "C": "nucl",
            "G": "nucl",
            "T": "nucl",
            "-": "gap",
        }
        default = "ambiguous"
    else:
        categories = {}
        for name in NAMES:
            categories[name] = "aa"
        categories["-"] = "gap"
        default = "ambiguous"

    categoryWidth = max(
        [len(category) for category in categories.values()] + [len(default)]
    )

    minLength = args.minLength
    concise = args.concise

    for index, read in enumerate(reads, start=1):
        counts: dict[str, int] = defaultdict(int)
        readLen = len(read)
        width = int(log10(readLen)) + 1
        if not concise:
            summary: list[str] = []
            append = summary.append
            offset = 1
        for category, count in sequenceCategoryLengths(
            read, categories, defaultCategory=default, minLength=minLength
        ):
            counts[category] += count
            if not concise:
                append(
                    "    %*d %-*s (offset %*d)"
                    % (width, count, categoryWidth, category, width, offset)
                )
                offset += count
        print("%d: %s (length %d)" % (index, read.id, readLen))
        for category in sorted(counts):
            count = counts[category]
            print(
                "  %-*s: %*d (%6.2f%%)"
                % (categoryWidth, category, width, count, count / readLen * 100.0)
            )
        if not concise:
            print("\n".join(summary))

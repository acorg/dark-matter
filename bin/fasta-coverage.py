#!/usr/bin/env python

import sys
from collections import defaultdict
import csv
from operator import itemgetter
from typing import Union

from dark.dna import AMBIGUOUS
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions

# These nucleotide codes will be listed first if --sortChars is used.
SPECIAL_CHARS = "ACGT-?N"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Print coverage information to standard output."
    )

    parser.add_argument(
        "--digits",
        default=4,
        type=int,
        help=(
            "The number of digits of precision to use for floating point "
            "outputs (coverage, gc fraction)."
        ),
    )

    parser.add_argument(
        "--noNames",
        dest="printNames",
        action="store_false",
        help="Do not print the ids of input sequences.",
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Treat all ambiguous nucleotides as not covered. Additional "
            "characters specified with --chars will also be counted as not "
            "covered."
        ),
    )

    parser.add_argument(
        "--gc", action="store_true", help="Include the GC percentage in the output."
    )

    parser.add_argument(
        "--reverse", "-r", action="store_true", help="Sort lines by decreasing value."
    )

    parser.add_argument(
        "--sortBy",
        help=(
            "What to sort lines by. Either one of id, coverage, length, gc, "
            "not-covered or else a single nucleotide character."
        ),
    )

    parser.add_argument(
        "--header",
        action="store_true",
        help="Write a header line and in cells print just the values.",
    )

    parser.add_argument("--sep", default="\t", help="The output separator character.")

    parser.add_argument(
        "--printCounts",
        action="store_true",
        help="Print the count for each character found in each read.",
    )

    parser.add_argument(
        "--sortChars", action="store_true", help="Sort the output columns."
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)
    coveredChars = set("AGCT") if args.strict else set(AMBIGUOUS)
    data = []
    sortBy = args.sortBy
    keepCounts = args.printCounts or args.gc or (sortBy and len(sortBy) == 1)

    for read in reads:
        notCoveredCount = 0
        if keepCounts:
            counts: dict[str, int] = defaultdict(int)
        for char in read.sequence.upper():
            notCoveredCount += char not in coveredChars
            if keepCounts:
                counts[char] += 1

        length = len(read)

        # Sanity check.
        if keepCounts:
            assert sum(counts.values()) == length

        coverage = round(
            (length - notCoveredCount) / length if length else 0.0, args.digits
        )

        readData = {
            "id": read.id,
            "length": len(read),
            "not-covered": notCoveredCount,
            "coverage": coverage,
        }

        if args.printCounts:
            readData["counts"] = counts

        if args.gc:
            gc = 0
            for char in "GC":
                if char in counts:
                    gc += counts[char]
            readData["gc"] = round(gc / length, args.digits)

        data.append(readData)

    if not data:
        sys.exit(0)

    writerow = csv.writer(sys.stdout, delimiter=args.sep).writerow

    if args.sortChars:
        allChars = set()
        for d in data:
            allChars.update(d["counts"])

        specials = [char for char in SPECIAL_CHARS if char in allChars]
        allCharsList = specials + sorted(allChars - set(specials))

    header = args.header

    if header:
        hline = ["Id"] if args.printNames else []
        hline.extend(("Length", "Not covered", "Coverage"))
        if args.gc:
            hline.append("GC")
        if args.printCounts and args.sortChars:
            hline.extend(allCharsList)
        writerow(hline)

    if sortBy:
        if sortBy in data[0]:
            key = itemgetter(sortBy)  # type: ignore
        else:
            if len(sortBy) > 1:
                print("Cannot sort by unknown column {sortBy!r}.", file=sys.stderr)
                sys.exit(1)

            sortBy = sortBy.upper()

            def key(d):  # type: ignore
                return d["counts"].get(sortBy, 0)

        data = sorted(data, key=key, reverse=args.reverse)

    for d in data:
        line: list[Union[str, int]] = [d["id"]] if args.printNames else []
        line.append(d["length"] if header else f"length:{d['length']}")
        line.append(d["not-covered"] if header else f"not-covered:{d['not-covered']}")
        line.append(d["coverage"] if header else f"coverage:{d['coverage']}")
        if args.gc:
            line.append(d["gc"] if header else f"gc:{d['gc']}")

        if args.printCounts:
            counts = d["counts"]
            if args.sortChars:
                for char in allCharsList:
                    if char in counts:
                        line.append(
                            counts[char] if header else f"{char}:{counts[char]}"
                        )
                    else:
                        line.append("")
            else:
                for char, count in sorted(counts.items()):
                    line.append(count if header else f"{char}:{count}")

        writerow(line)

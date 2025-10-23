#!/usr/bin/env python

import argparse
import sys

from dark.reads import DNARead, addFASTACommandLineOptions, parseFASTACommandLineOptions

FORWARD, REVERSE, COMPLEMENT, REVERSE_COMPLEMENT = 0, 1, 2, 3

SUFFIX = {
    FORWARD: "",
    REVERSE: " (r)",
    COMPLEMENT: " (c)",
    REVERSE_COMPLEMENT: " (rc)",
}


def getParser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin and a sequence (or genome region) to find, write "
            "information about (potentially overlapping) matches (and, optionally, "
            "reversed, complemented, or reversed complement matches), to standard "
            "output. Offsets are inclusive and 1-based unless --zeroBased is used "
            "(in which case they are 0-based closed/open offsets as used in Python). "
        )
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("--pattern", "--sequence", help="The sequence to look for.")

    group.add_argument(
        "--region",
        help=(
            "A region of the form A-B, where A and B are numeric sites, with "
            "A <= B. A and B will be treated as 1-based offsets unless "
            "--zeroBased is used. This can be used to find repeat copies of an "
            "input sequence region. I.e., for each input sequence, the specified "
            "region of the sequence will be searched for in the rest of the sequence."
        ),
    )

    parser.add_argument(
        "--reverseComplement",
        "--rc",
        action="store_true",
        help=(
            "Also search for the reverse complement. Reverse complement matches "
            f"are indicated with an{SUFFIX[REVERSE_COMPLEMENT]} suffix."
        ),
    )

    parser.add_argument(
        "--reverse",
        action="store_true",
        help=(
            f"Also search for the reversed sequence. Reverse matches are "
            f"indicated with an{SUFFIX[REVERSE]} suffix."
        ),
    )

    parser.add_argument(
        "--complement",
        action="store_true",
        help=(
            f"Also search for the complemented sequence. Complement matches "
            f"are indicated with a{SUFFIX[COMPLEMENT]} suffix."
        ),
    )

    parser.add_argument(
        "--end",
        action="store_true",
        help="Print the index of the end of the match, not the start.",
    )

    parser.add_argument(
        "--offsetsOnly",
        action="store_true",
        help=(
            "Only print the matching offset (or offsets, separated by the --sep value)."
        ),
    )

    parser.add_argument(
        "--checkMatchCount",
        type=int,
        help=(
            "The expected match count. This can be used to force an error if the "
            "number of occurrences of the search pattern found is not as expected."
        ),
    )

    parser.add_argument(
        "--caseSensitive",
        action="store_true",
        help="Do not ignore upper/lower case when matching.",
    )

    parser.add_argument(
        "--printSequences",
        action="store_true",
        help=(
            f"Print the sequences searched for. Reversed, complemented, and "
            f"reverse-complemented sequences are indicated by "
            f"{SUFFIX[REVERSE]},{SUFFIX[COMPLEMENT]}, "
            f"and{SUFFIX[REVERSE_COMPLEMENT]}, respectively."
        ),
    )

    parser.add_argument(
        "--printTrivialMatch",
        action="store_true",
        help=(
            "If --region is used, there is a guaranteed match at that "
            "location. Do not include that match in the printed matches "
            "unless this option is given."
        ),
    )

    parser.add_argument(
        "--sep", default="\t", help="The output field separator string (default=TAB)."
    )

    parser.add_argument(
        "--zeroBased",
        action="store_true",
        help=(
            "Return Python-style zero-based closed/open match offsets. When --end is "
            "used, the output offset will be for the character that follows the end of "
            "the match (i.e., it is not part of the match)."
        ),
    )

    parser.add_argument(
        "--gapCharacter",
        default="-",
        help="The gap character (only used if --ignoreGaps is used).",
    )

    parser.add_argument(
        "--ignoreGaps",
        action="store_true",
        help=(
            "Ignore gaps when matching. This is useful for finding a subsequence of "
            "an aligned sequence. Note that the reported index into input sequences "
            "will be for the (possibly gapped) original input sequence."
        ),
    )

    addFASTACommandLineOptions(parser)

    return parser


def parseRegion(region: str, zeroBased: bool) -> tuple[int, int]:
    """
    Parse a region string like 34-79 and return the two numbers.
    """
    lo, hi = map(int, region.split("-"))

    if zeroBased:
        if not (0 <= lo <= hi):
            sys.exit(
                f"Invalid range {region!r}. Give two non-negative numbers "
                "in non-decreasing order.",
            )
    else:
        if not (0 < lo <= hi):
            extra = (
                " Or use --zeroBased if you really want to specify a zero offset."
                if lo == 0
                else ""
            )

            sys.exit(
                f"Invalid range {region!r}. Give two positive numbers in "
                f"non-decreasing order.{extra}",
            )

        # Use zero-based indices.
        lo -= 1

    return lo, hi


def getPatterns(
    pattern: str, reverse: bool, complement: bool, reverseComplement: bool
) -> list[tuple[str, int]]:
    """
    Get a list of patterns (and their orientations) to search for.
    """
    # We always look for the given pattern in the forward direction.
    result = [(pattern, FORWARD)]

    if reverse:
        result.append((pattern[::-1], REVERSE))

    rc = DNARead("id", pattern).reverseComplement().sequence

    if reverseComplement:
        result.append((rc, REVERSE_COMPLEMENT))

    if complement:
        # The complement is just the reverse of the reverse complement.
        result.append((rc[::-1], COMPLEMENT))

    return result


def main() -> None:
    args = getParser().parse_args()

    if args.region:
        lo, hi = parseRegion(args.region, args.zeroBased)
        patterns = []  # Will be set below, extracted from each target read.
    else:
        lo = hi = -1  # Unused.
        patterns = getPatterns(
            args.pattern, args.reverse, args.complement, args.reverseComplement
        )

    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        if args.region:
            patterns = getPatterns(
                read.sequence[lo:hi],
                args.reverse,
                args.complement,
                args.reverseComplement,
            )

        matches: set[tuple[int, int]] = set()

        for pattern, orientation in patterns:
            # Look for as many occurrences of this pattern as can be found in the read
            # sequence.
            start = 0
            while True:
                index = read.find(
                    pattern,
                    start=start,
                    end=args.end,
                    caseSensitive=args.caseSensitive,
                    ignoreGaps=args.ignoreGaps,
                    gapCharacter=args.gapCharacter,
                )
                if index == -1:
                    # No more matches of this pattern.
                    break
                else:
                    if args.region is None or args.printTrivialMatch or index != lo:
                        assert (index, orientation) not in matches
                        matches.add((index, orientation))
                    start = index + 1

        if args.checkMatchCount is not None and len(matches) != args.checkMatchCount:
            sys.exit(
                f"{len(matches)} {'match was' if len(matches) == 1 else 'matches were'} "  # noqa: E501
                f"found for sequence {read.id!r}, but {args.checkMatchCount} "
                f"{'match was' if args.checkMatchCount == 1 else 'matches were'} "
                "expected."
            )

        if matches:
            adjust = 0 if args.end or args.zeroBased else 1

            if args.offsetsOnly:
                fields = [str(index + adjust) for index, _ in sorted(matches)]
            else:
                fields = [
                    read.id,
                    f"{len(matches)} match{'' if len(matches) == 1 else 'es'}",
                    ", ".join(
                        str(index + adjust) + SUFFIX[orientation]
                        for index, orientation in sorted(matches)
                    ),
                ]

                if args.printSequences:
                    for pattern, orientation in patterns:
                        fields.append(f"sequence{SUFFIX[orientation]} = {pattern}")

            print(args.sep.join(fields))


if __name__ == "__main__":
    main()

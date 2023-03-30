#!/usr/bin/env python

import sys

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions, DNARead)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin and a sequence (or genome region) "
            "to find, write information about matches (and, "
            "optionally, reverse complement matches), to standard "
            "ouput."
        )
    )

    parser.add_argument(
        "--reverseComplement",
        "--rc",
        action="store_true",
        help="Also search for the reverse complement.",
    )

    parser.add_argument(
        "--caseSensitive",
        action="store_true",
        help="Do not ignore upper/lower case when matching.",
    )

    parser.add_argument(
        "--printSequence", action="store_true",
        help="Print the sequence searched for."
    )

    parser.add_argument(
        "--includeShortGenomes",
        action="store_true",
        help=(
            "Try to match region sequences that are shorter than "
            "the length of the given region (assuming they are "
            "non-empty) due to the input sequence being too short."
        ),
    )

    parser.add_argument(
        "--noShortGenomeWarning",
        dest="reportShortGenomes",
        action="store_false",
        help=(
            "Do not report (to standard error) when sequences are shorter "
            "than the length of the given region due to being too short."
        ),
    )

    parser.add_argument(
        "--printReverseComplement",
        "--printRc",
        action="store_true",
        help=(
            "Print the reverse complemented region searched for (only "
            "valid with --reverseComplement)."
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
        "--sep", default="\t", help="The field separator string (default=TAB)."
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "--region",
        help=(
            "A 1-based region of the form A-B, where A and B are numeric "
            "sites, with A <= B. Both indices are included."
        ),
    )

    group.add_argument("--sequence", help="The sequence to look for.")

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()

    if args.printReverseComplement and not args.reverseComplement:
        print(
            "The --printReverseComplement option only makes sense if you "
            "also use --reverseComplement.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.region:
        lo, hi = map(int, args.region.split("-"))

        if not (0 < lo <= hi):
            print(
                f"Invalid range {args.region!r}. Give two non-negative "
                f"numbers in non-decreasing order.",
                file=sys.stderr,
            )
            sys.exit(1)

        # Use zero-based indices.
        lo -= 1
    else:
        target = args.sequence

    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        if args.region:
            target = read.sequence[lo:hi]
            if (len(target) < hi - lo and not args.includeShortGenomes) or len(
                target
            ) == 0:
                if args.reportShortGenomes:
                    print(
                        f"The {read.id!r} sequence is too short ({len(read)}) "
                        f"for a full match of a region of length {hi - lo} "
                        f"starting at site {lo + 1} in the genome.",
                        file=sys.stderr,
                    )
                continue

        if args.reverseComplement:
            targets = (target,
                       DNARead("id", target).reverseComplement().sequence)
            rcs = False, True
        else:
            targets = (target,)
            rcs = (False,)

        sequence = read.sequence

        if not args.caseSensitive:
            # list is needed in the following as we may use 'targets'
            # twice.
            targets = list(map(str.upper, targets))

        matches = set()

        for thisTarget, rc in zip(targets, rcs):
            start = 0
            while True:
                index = sequence.find(thisTarget, start)
                if index == -1:
                    break
                else:
                    if (args.region is None or args.printTrivialMatch or
                            index != lo):
                        assert (index, rc) not in matches
                        matches.add((index, rc))
                    start = index + 1

        if matches:
            fields = [
                read.id,
                "%d match%s" % (len(matches),
                                "" if len(matches) == 1 else "es"),
                ", ".join(
                    str(index + 1) + (" (rc)" if rc else "")
                    for index, rc in sorted(matches)
                ),
            ]

            for thisTarget, rc in zip(targets, rcs):
                if rc:
                    if args.printReverseComplement:
                        fields.append("sequence (rc) = " + thisTarget)
                else:
                    if args.printSequence:
                        fields.append("sequence = " + thisTarget)

            print(args.sep.join(fields))

#!/usr/bin/env python

import sys
import argparse

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions, DNARead


FORWARD, REVERSE, COMPLEMENT, REVERSE_COMPLEMENT = 0, 1, 2, 3

SUFFIX = {FORWARD: "", REVERSE: " (r)", COMPLEMENT: " (c)", REVERSE_COMPLEMENT: " (rc)"}


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
    help="Also search for the reverse complement. Reverse complement matches "
    f"are indicated with an{SUFFIX[REVERSE_COMPLEMENT]} suffix.",
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

parser.add_argument(
    "--zeroBased", action="store_true", help="Print zero-based matched offsets."
)

group = parser.add_mutually_exclusive_group(required=True)

group.add_argument(
    "--region",
    help=(
        "A region of the form A-B, where A and B are numeric sites, with "
        "A <= B. A and B will be treated as 1-based offsets unless "
        "--zeroBased is used."
    ),
)

group.add_argument("--sequence", help="The sequence to look for.")

addFASTACommandLineOptions(parser)
args = parser.parse_args()

if args.region:
    lo, hi = map(int, args.region.split("-"))

    if args.zeroBased:
        if not (0 <= lo <= hi):
            print(
                f"Invalid range {args.region!r}. Give two non-negative "
                f"numbers in non-decreasing order.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        if not (0 < lo <= hi):
            extra = (
                " Or use --zeroBased if you really want to specify " "a zero offset."
                if lo == 0
                else ""
            )

            print(
                f"Invalid range {args.region!r}. Give two positive "
                f"numbers in non-decreasing order.{extra}",
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
        if (len(target) < hi - lo and not args.includeShortGenomes) or len(target) == 0:
            if args.reportShortGenomes:
                print(
                    f"The {read.id!r} sequence is too short ({len(read)}) "
                    f"for a full match of a region of length {hi - lo} "
                    f"starting at site {lo + 1} in the genome.",
                    file=sys.stderr,
                )
            continue

    # Always look for the given target.
    targets = [target]
    orientations = [FORWARD]

    # And look for target variants according to Boolean command-line options.
    if args.reverse:
        targets.append(target[::-1])
        orientations.append(REVERSE)

    if args.complement:
        targets.append(DNARead("id", target).reverseComplement().sequence[::-1])
        orientations.append(COMPLEMENT)

    if args.reverseComplement:
        targets.append(DNARead("id", target).reverseComplement().sequence)
        orientations.append(REVERSE_COMPLEMENT)

    sequence = read.sequence

    if not args.caseSensitive:
        # list is needed in the following as we may use 'targets'
        # twice.
        targets = list(map(str.upper, targets))

    matches = set()

    for thisTarget, orientation in zip(targets, orientations):
        start = 0
        while True:
            index = sequence.find(thisTarget, start)
            if index == -1:
                break
            else:
                if args.region is None or args.printTrivialMatch or index != lo:
                    assert (index, orientation) not in matches
                    matches.add((index, orientation))
                start = index + 1

    if matches:
        adjust = 0 if args.zeroBased else 1
        fields = [
            read.id,
            "%d match%s" % (len(matches), "" if len(matches) == 1 else "es"),
            ", ".join(
                str(index + adjust) + SUFFIX[orientation]
                for index, orientation in sorted(matches)
            ),
        ]

        if args.printSequences:
            for thisTarget, orientation in zip(targets, orientations):
                fields.append(("sequence%s = " % SUFFIX[orientation]) + thisTarget)

        print(args.sep.join(fields))

#!/usr/bin/env python

import sys
from typing import Optional

from dark.errors import ReadLengthsNotIdenticalError
from dark.reads import (
    addFASTACommandLineOptions,
    parseFASTACommandLineOptions,
    Reads,
    Read,
)


def printHeader(variableSites, args, baseOffset):
    n = len(variableSites)
    print(
        "%d site%s %svariable (threshold for homogeneity: %.3f)."
        % (
            n,
            " was" if n == 1 else "s were",
            "confirmed " if args.confirm else "",
            args.homogeneous,
        ),
        file=sys.stderr,
    )
    if args.reference:
        print(
            "Offsets adjusted by %d, relative to sequence %r."
            % (baseOffset, args.reference),
            file=sys.stderr,
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on standard input, write reads with only the "
            "variable sites to standard output. If no sites are variable, "
            "nothing is written to standard output."
        )
    )

    parser.add_argument(
        "--printSites",
        action="store_true",
        default=False,
        help=(
            "Print the variable site locations to standard error along "
            "with the nucleotides found at the site, the homogeneity "
            "fraction of the site and the reference base (if a reference "
            "is given)."
        ),
    )

    parser.add_argument(
        "--sitesOnly",
        action="store_true",
        default=False,
        help=(
            "Only print the (comma separated) variable sites to standard "
            "output. This output can then be used directly as an argument "
            "for --keepSites to filter-fasta.py"
        ),
    )

    parser.add_argument(
        "--homogeneous",
        default=1.0,
        type=float,
        help=(
            "If the most-common nucleotide frequency at a site is at least "
            "this value, the site will be considered homogeneous."
        ),
    )

    parser.add_argument(
        "--confirm",
        action="store_true",
        default=False,
        help=(
            "Only keep sites where there is confirm variation (i.e., "
            "ambiguous sites that are compatible with there being no "
            "variation are not included)."
        ),
    )

    parser.add_argument(
        "--reference",
        help=(
            "A sequence id. Site offsets will be reported relative to this "
            "sequence. That means the number of gap (hyphen) chars at the "
            "start of this sequence will be counted and offsets will be "
            "incremented by that amount (and lower offsets in sequences "
            "that start before the specified sequence will be ignored)."
        ),
    )

    parser.add_argument(
        "--unknownAreAmbiguous",
        action="store_true",
        default=False,
        help=(
            "Any unknown character (e.g., a '-' gap or '?' unknown base) "
            "will be treated as being fully ambiguous (i.e., could be any "
            "of ACGT). Otherwise, all unknown characters are counted "
            "as '-' characters."
        ),
    )

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = Reads(list(parseFASTACommandLineOptions(args)))

    if not reads:
        sys.exit(0)

    reference: Optional[Read]

    if args.reference:
        for read in reads:
            if read.id == args.reference:
                break
        else:
            print(
                "Could not find --reference sequence %r." % args.reference,
                file=sys.stderr,
            )
            sys.exit(1)
        baseOffset = len(read.sequence) - len(read.sequence.lstrip("-"))
        reference = read
    else:
        baseOffset = 0
        reference = None

    try:
        variableSites = reads.variableSites(
            confirm=args.confirm,
            homogeneityLevel=args.homogeneous,
            unknownAreAmbiguous=args.unknownAreAmbiguous,
        )
    except ReadLengthsNotIdenticalError:
        print("Input sequences are not all the same length!", file=sys.stderr)
        sys.exit(2)

    if variableSites:
        toDelete = set()
        if args.printSites:
            for site, counts in variableSites.items():
                if site >= baseOffset:
                    ref = (" (ref %s)" % reference.sequence[site]) if reference else ""
                    print(
                        "%d: %s%s" % (site + 1 - baseOffset, counts, ref),
                        file=sys.stderr,
                    )
                else:
                    toDelete.add(site)

        for site in toDelete:
            del variableSites[site]

    if variableSites:
        if args.sitesOnly:
            print(
                ",".join(
                    map(lambda site: str(site + 1 - baseOffset), sorted(variableSites))
                )
            )
        else:
            saveAs = "fasta" if args.fasta else "fastq"
            reads.filter(keepSites=set(variableSites)).save(sys.stdout, saveAs)
            printHeader(variableSites, args, baseOffset)
    else:
        print(
            "No sites were %svariable (threshold for homogeneity: %.3f)."
            % ("confirmed " if args.confirm else "", args.homogeneous),
            file=sys.stderr,
        )

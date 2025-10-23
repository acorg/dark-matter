#!/usr/bin/env python

import argparse
import multiprocessing
import sys
from math import log10

from dark.aligners import MAFFT_DEFAULT_ARGS, NEEDLE_DEFAULT_ARGS, align
from dark.dna import AMBIGUOUS, compareDNAReads, matchToString
from dark.reads import Reads, addFASTACommandLineOptions, parseFASTACommandLineOptions
from dark.utils import parseRangeExpression

MAFFT_ALGORITHMS_URL = (
    "https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html"
)


def getArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compare two DNA sequences.",
    )

    parser.add_argument(
        "--index1",
        type=int,
        default=1,
        help="The (1-based) index in the input of the first sequence.",
    )

    parser.add_argument(
        "--index2",
        type=int,
        default=2,
        help="The (1-based) index in the input of the second sequence.",
    )

    parser.add_argument(
        "--align",
        action="store_true",
        help=(
            "If given, use mafft (the default) or needle (according to the "
            "algorithm selected by --aligner) to align the two sequences."
        ),
    )

    parser.add_argument(
        "--aligner",
        default="mafft",
        choices=("edlib", "mafft", "needle"),
        help="The alignment algorithm to use.",
    )

    parser.add_argument(
        "--alignerOptions",
        help=(
            "Optional arguments to pass to the alignment algorithm. If the "
            f"aligner is 'mafft', the default options are {MAFFT_DEFAULT_ARGS!r}. "
            f"If 'needle', the default is {NEEDLE_DEFAULT_ARGS!r}. Do not try to set "
            "the number of threads here; use the --threads argument instead. If you "
            f"are using mafft, see {MAFFT_ALGORITHMS_URL} for some possible option "
            "combinations."
        ),
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help=(
            "The number of threads to use when running the aligner (if --align "
            "is used and the alignment algorithm can make use of multiple "
            "threads (mafft can, needle cannot))."
        ),
    )

    parser.add_argument(
        "--alignmentFile", help="The file to save the alignment to (implies --align)."
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="If given, do not allow ambiguous nucleotide symbols to match.",
    )

    parser.add_argument(
        "--quiet",
        dest="verbose",
        action="store_false",
        help="Do not print information about aligning.",
    )

    parser.add_argument(
        "--noGapLocations",
        dest="includeGapLocations",
        action="store_false",
        help="Do not indicate the (1-based) locations of sequence gaps.",
    )

    parser.add_argument(
        "--noNoCoverageLocations",
        dest="includeCoverageLocations",
        action="store_false",
        help="Do not indicate the (1-based) locations of no coverage.",
    )

    parser.add_argument(
        "--sites",
        help=(
            "Specify (1-based) sequence sites to keep. All other sites will "
            "be ignored. The sites must be given in the form e.g., 24,100-200,260."
        ),
    )

    parser.add_argument(
        "--showDiffs",
        action="store_true",
        help="Print (1-based) sites where the sequence nucleotides differ.",
    )

    parser.add_argument(
        "--showAmbiguous",
        action="store_true",
        help=(
            "Print (1-based) sites where either sequence has an ambiguous "
            "nucleotide code."
        ),
    )

    parser.add_argument(
        "--includeAmbiguousMatches",
        action="store_true",
        help=(
            "Print (1-based) sites of ambiguous matches. The output gives the "
            "site, the base from the first sequence, then the base from the "
            "second sequence."
        ),
    )

    parser.add_argument(
        "--includeNonGapMismatches",
        action="store_true",
        help=(
            "Print (1-based) sites of mismatches that do not involve a gap. The "
            "output gives the site, the base from the first sequence, then the "
            "base from the second sequence."
        ),
    )

    parser.add_argument(
        "--includeDifferenceCounts",
        action="store_true",
        help=(
            "Include the counts (i.e., number of occurrences) of all differences. "
            "Note that mismatches involving areas of no coverage or where one read "
            "is shorter than another do not count as differences."
        ),
    )

    parser.add_argument(
        "--includeDifferenceLocations",
        action="store_true",
        help=(
            "Include the locations (and counts) (i.e., number of occurrences) "
            "of all differences. Note that mismatches involving areas of no coverage "
            "or where one read is shorter than another do not count as differences."
        ),
    )

    parser.add_argument(
        "--gapChars",
        default="-",
        metavar="CHARS",
        help=(
            "The sequence characters that should be considered to be gaps. "
            "These characters will be ignored in computing sequence lengths "
            "and identity fractions."
        ),
    )

    parser.add_argument(
        "--noCoverageChars",
        metavar="CHARS",
        help=(
            "The sequence characters that indicate lack of coverage. "
            "These characters will be ignored in identity fractions."
        ),
    )

    addFASTACommandLineOptions(parser)

    return parser.parse_args()


def main() -> None:
    args = getArgs()

    keepSequences = set([args.index1 - 1, args.index2 - 1])

    reads = list(parseFASTACommandLineOptions(args).filter(keepSequences=keepSequences))

    if len(reads) == 1:
        if len(keepSequences) == 1:
            # This is odd but not necessarily an error. For some reason, they're asking
            # for a comparison of a sequence with itself.
            if args.verbose:
                print(
                    "Warning: you are comparing a sequence with itself.",
                    file=sys.stderr,
                )

            reads = Reads([reads[0], reads[0]])
        else:
            sys.exit("Could not find both requested sequence indices. Exiting.")
    elif len(reads) != 2:
        sys.exit("Could not find both requested sequence indices. Exiting.")

    if args.alignmentFile:
        args.align = True

    if args.align:
        len1, len2 = map(len, reads)
        if len1 == len2:
            print("Pre-alignment, sequence lengths were identical: %d" % len1)
        else:
            print(
                "Pre-alignment, sequence lengths: %d, %d (difference %d)"
                % (len1, len2, abs(len1 - len2))
            )

        print("  Gaps:")
        print("    Id: %s %d" % (reads[0].id, reads[0].sequence.count("-")))
        print("    Id: %s %d" % (reads[1].id, reads[1].sequence.count("-")))

        reads = align(
            reads, args.aligner, args.alignerOptions, args.verbose, args.threads
        )

        if args.alignmentFile:
            assert reads.save(args.alignmentFile) == 2

    offsets = (
        parseRangeExpression(args.sites, convertToZeroBased=True)
        if args.sites
        else None
    )

    read1, read2 = reads
    len1, len2 = map(len, reads)
    identicalLengths = len1 == len2

    # Sanity check.
    if args.align:
        assert identicalLengths

    match = compareDNAReads(
        read1,
        read2,
        matchAmbiguous=(not args.strict),
        offsets=offsets,
        gapChars=args.gapChars,
        noCoverageChars=args.noCoverageChars,
    )

    x = "Post-alignment, sequence" if args.align else "Sequence"
    if identicalLengths:
        print("%s lengths are identical: %s" % (x, len1))
    else:
        print("%s lengths: %d, %d (difference %d)" % (x, len1, len2, abs(len1 - len2)))

    print(
        matchToString(
            match,
            read1,
            read2,
            matchAmbiguous=(not args.strict),
            offsets=offsets,
            includeGapLocations=args.includeGapLocations,
            includeNoCoverageLocations=args.includeCoverageLocations,
            includeAmbiguousMatches=args.includeAmbiguousMatches,
            includeNonGapMismatches=args.includeNonGapMismatches,
            includeDifferenceLocations=args.includeDifferenceLocations,
            includeDifferenceCounts=args.includeDifferenceCounts,
        )
    )

    if args.showDiffs:
        # Print all sites where the sequences differ.
        width = int(log10(max(len1, len2))) + 1
        headerPrinted = False
        for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence), start=1):
            if a != b:
                if not headerPrinted:
                    print("Differences (site, %s, %s):" % (read1.id, read2.id))
                    headerPrinted = True
                print("  %*d %s %s" % (width, site, a, b))

        if not headerPrinted:
            print("No sequence differences found.")

    if args.showAmbiguous:
        width = int(log10(max(len1, len2))) + 1
        headerPrinted = False
        for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence), start=1):
            if len(AMBIGUOUS.get(a, "")) > 1 or len(AMBIGUOUS.get(b, "")) > 1:
                if not headerPrinted:
                    print("Ambiguities (site, %s, %s):" % (read1.id, read2.id))
                    headerPrinted = True
                print("  %*d %s %s" % (width, site, a, b))

        if not headerPrinted:
            print("No sequence ambiguities found.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python

"""
Read DNA or AA FASTA from stdin and print AA ORF-only FASTA to stdout. I.e.,
the output FASTA will contain one sequence for each ORF found in each
sequence on stdin.

If a minimal ORF length is given, only print ORFs of at least that length
unless we are asked to allow ORFs that are open at at least one end.

Note that no start or stop codons will appear in the output. If you are
wanting to simply translate DNA FASTA to AA FASTA on a sequence-by-sequence
basis leaving in start/stop codons, use dna-to-aa.py instead.
"""

import sys
from os.path import basename
import argparse

from Bio.Data.CodonTable import TranslationError

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert DNA or AA FASTA to AA ORF-only FASTA",
        epilog=(
            "Given DNA or AA FASTA on stdin, output AA ORF-only FASTA to "
            "stdout. Optionally, filter by minimum required ORF length "
            "and allow the output of short ORFs that are open."
        ),
    )

    parser.add_argument(
        "--allowOpenORFs",
        default=False,
        action="store_true",
        help=(
            "If True, ORFs that do not meet the length requirement "
            "(as specified by --minORFLength) will be output as long as "
            "they are open."
        ),
    )

    parser.add_argument(
        "--minORFLength",
        metavar="LEN",
        type=int,
        default=None,
        help="Only ORFs of at least this length will be written to stdout.",
    )

    parser.add_argument(
        "--maxORFLength",
        metavar="LEN",
        type=int,
        default=None,
        help=("Only ORFs of a maximum of this length will be written to " "stdout."),
    )

    parser.add_argument(
        "--kozakOnly",
        default=False,
        action="store_true",
        help=(
            "Only ORFs that also have a Kozak consensus will be written to "
            "stdout. Only applicable if DNA reads are given."
        ),
    )

    parser.add_argument(
        "--kozakInfoFile",
        help=(
            "Filename to which all Kozak consensus information is written. "
            "Only applicable if DNA reads are given."
        ),
    )

    addFASTACommandLineOptions(parser)

    args = parser.parse_args()
    allowOpenORFs = args.allowOpenORFs
    minORFLength = args.minORFLength
    maxORFLength = args.maxORFLength
    kozakInfoFile = args.kozakInfoFile
    kozakOnly = args.kozakOnly
    reads = parseFASTACommandLineOptions(args)
    aa = args.readClass in (
        "AARead",
        "AAReadORF",
        "AAReadWithX",
        "SSAARead",
        "SSAAReadWithX",
        "TranslatedRead",
    )

    if aa:

        def translations(read):
            return (read,)

    else:

        def translations(read):
            return read.translations()

    if kozakInfoFile or kozakOnly:
        if aa:
            print(
                "Kozak sequences cannot be computed from aa sequences.", file=sys.stderr
            )
            sys.exit(1)
        else:
            from dark.dna import findKozakConsensus

        def writeToKozakOut(kozakread, kozakfp):
            """
            Writes out information about a Kozak sequence stored in kozakread
            to the file name given in kozakfp.
            """
            print("Read ID: " + kozakread.id, file=kozakfp)
            print("Kozak sequence: " + kozakread.sequence, file=kozakfp)
            print(
                "Offset of the Kozak sequence A: " + str(kozakread.start + 6),
                file=kozakfp,
            )
            print(
                "Quality of the Kozak consensus: " + str(kozakread.kozakQuality) + " %",
                file=kozakfp,
            )

    def findORFs(kozakfp):
        for read in reads:
            try:
                for translation in translations(read):
                    for orf in translation.ORFs(allowOpenORFs):
                        # Check the length requirements, if any.
                        if (minORFLength is None or len(orf) >= minORFLength) and (
                            maxORFLength is None or len(orf) <= maxORFLength
                        ):
                            if kozakOnly:
                                for kozakread in findKozakConsensus(read):
                                    start = orf.start * 3
                                    if start + 1 <= kozakread.stop <= start + 3:
                                        print(orf.toString("fasta"), end="")
                            else:
                                print(orf.toString("fasta"), end="")
                            if kozakfp:
                                for kozakread in findKozakConsensus(read):
                                    start = orf.start * 3
                                    if start + 1 <= kozakread.stop <= start + 3:
                                        writeToKozakOut(kozakread, kozakfp)

            except TranslationError as error:
                print(
                    "Could not translate read %r sequence %r (%s)."
                    % (read.id, read.sequence, error),
                    file=sys.stderr,
                )
                if not aa:
                    print(
                        "Did you forget to run "
                        '%s with "--readClass AARead"?' % (basename(sys.argv[0])),
                        file=sys.stderr,
                    )
                sys.exit(1)

    if kozakInfoFile:
        with open(kozakInfoFile, "a+") as kozakfp:
            findORFs(kozakfp)
    else:
        findORFs(None)

#!/usr/bin/env python

import sys
import argparse
from itertools import chain

from dark.civ.proteins import SqliteIndex
from dark.errors import NoSuchGenomeError
from dark.filter import (
    addFASTAFilteringCommandLineOptions,
    parseFASTAFilteringCommandLineOptions,
)
from dark.genbank import GenomeRanges
from dark.genomes import GenomeProteinInfo
from dark.reads import Reads
from dark.sam import SAMFilter, samReferences
from dark.utils import pct


def summarize(gpi, sortOn, minReadOffsetCount):
    """
    Print a summary of the genome proteins.

    @param gpi: A C{GenomeProteinInfo} instance.
    @param sortOn: How to sort proteins for output. One of 'coverage',
        'depth', 'name', 'offset', or 'readCount'.
    @param minReadOffsetCount: The minimum number of reads offsets that must
        overlap a protein for the read to be considered as sufficiently
        intersecting the protein.
    """
    genome = gpi.genome

    print("Summary of %s (%s):" % (genome["name"], genome["accession"]))
    print("  Length: %d" % genome["length"])
    print("  Protein count: %d" % genome["proteinCount"])
    print("  Total protein offsets: %s" % (pct(len(gpi.offsets), genome["length"])))

    if gpi.samFiles:
        print("  SAM files analyzed: %d" % len(gpi.samFiles))
        for i, filename in enumerate(gpi.samFiles, start=1):
            print("    %d: %s" % (i, filename))
    else:
        return

    print("  Whole genome coverage (not just proteins):")
    print("    Reads matching genome: %d" % len(gpi.readIdsMatchingGenome))
    print(
        "    Covered genome offsets: %s"
        % (pct(len(gpi.coveredOffsetCount), genome["length"]))
    )
    print(
        "    Average depth across genome: %.3f"
        % (sum(gpi.coveredOffsetCount.values()) / genome["length"])
    )

    coveredProteinOffsetCount = coveredProteinBasesCount = 0
    for offset in gpi.offsets:
        coveredProteinOffsetCount += bool(gpi.coveredOffsetCount[offset])
        coveredProteinBasesCount += gpi.coveredOffsetCount[offset]

    print("  Total protein coverage (irrespective of minReadOffsetCount):")
    print("    Reads matching proteins: %d" % len(gpi.readIdsForAllProteins()))
    print(
        "    Proteins with any coverage: %s"
        % pct(len(gpi.coveredProteins), genome["proteinCount"])
    )
    print(
        "    Covered protein offsets: %s"
        % (pct(coveredProteinOffsetCount, len(gpi.offsets)))
    )
    print(
        "    Average depth across proteins: %.3f"
        % (coveredProteinBasesCount / len(gpi.offsets))
    )

    if sortOn == "name":

        def key(proteinAccession):
            return gpi.proteins[proteinAccession]["name"]

        reverse = False
    elif sortOn == "offset":

        def key(proteinAccession):
            return GenomeRanges(gpi.proteins[proteinAccession]["offsets"]).ranges[0][0]

        reverse = False
    elif sortOn == "readCount":

        def key(proteinAccession):
            coverage = gpi.proteinCoverageInfo(proteinAccession)
            return len(coverage["readIds"])

        reverse = True
    elif sortOn == "coverage":

        def key(proteinAccession):
            coverage = gpi.proteinCoverageInfo(proteinAccession)
            return coverage["coveredOffsets"] / coverage["ntLength"]

        reverse = True
    elif sortOn == "depth":

        def key(proteinAccession):
            coverage = gpi.proteinCoverageInfo(proteinAccession)
            return coverage["totalBases"] / coverage["ntLength"]

        reverse = True

    if minReadOffsetCount is None:
        print("  Proteins covered (no minReadOffsetCount):")
    else:
        print("  Proteins covered (minReadOffsetCount=%d):" % minReadOffsetCount)

    proteinCount = 0
    for proteinAccession in sorted(gpi.coveredProteins, key=key, reverse=reverse):
        protein = gpi.proteins[proteinAccession]

        coverage = gpi.proteinCoverageInfo(proteinAccession, minReadOffsetCount)

        readCount = len(coverage["readIds"])

        if readCount:
            proteinCount += 1

            print(
                "    %d: %s (%d AA, %d nt with stop codon, %s)"
                % (
                    proteinCount,
                    protein["product"],
                    protein["length"],
                    protein["length"] * 3 + 3,
                    protein["accession"],
                )
            )

            print("      Read count: %d" % readCount)

            print(
                "      Covered offsets: %s"
                % (pct(coverage["coveredOffsets"], coverage["ntLength"]))
            )

            print(
                "      Average depth: %.3f"
                % (coverage["totalBases"] / coverage["ntLength"])
            )

            print("      Offsets: %s" % protein["offsets"])

    print(
        "  Proteins matched: %s (sorted by %s):"
        % (pct(proteinCount, genome["proteinCount"]), sortOn)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Print SAM/BAM file protein match statistics.",
    )

    parser.add_argument(
        "--proteinGenomeDatabase",
        required=True,
        help=(
            "The filename of an Sqlite3 database holding protein and "
            "genome information, as built by make-protein-database.py"
        ),
    )

    parser.add_argument(
        "--progress",
        default=False,
        action="store_true",
        help="Print progress info to standard error.",
    )

    parser.add_argument(
        "--sortOn",
        default="readCount",
        choices=("coverage", "depth", "name", "offset", "readCount"),
        help="How to sort proteins for output.",
    )

    parser.add_argument(
        "--minReadOffsetCount",
        type=int,
        help=(
            "The minimum number of reads offsets that must overlap a "
            "protein for the read to be considered as sufficiently "
            "intersecting the protein. Use this to prevent reads that "
            "just overlap the protein in a very small number offsets "
            "from being counted."
        ),
    )

    parser.add_argument(
        "--skipTranslationChecks",
        dest="checkTranslations",
        action="store_false",
        default=True,
        help=(
            "Skip the sanity check that database protein sequences can all "
            "be translated from the database genome sequence."
        ),
    )

    addFASTAFilteringCommandLineOptions(parser)
    SAMFilter.addFilteringOptions(
        parser,
        samfileIsPositional=False,
        samfileAction="append",
        samfileNargs="*",
        samfileRequired=False,
    )

    args = parser.parse_args()

    samfiles = list(chain.from_iterable(args.samfile)) if args.samfile else []

    if samfiles:
        if args.referenceId:
            referenceIds = args.referenceId
        else:
            # If all SAM files have just one reference and they're all the
            # same, use that. Else complain.
            referenceIds = set()
            for filename in samfiles:
                referenceIds.update(samReferences(filename))

            if len(referenceIds) != 1:
                print(
                    "No reference id(s) specified with --referenceId, and "
                    "the given SAM/BAM files do not contain exactly one "
                    "(identical) reference. Please use --referenceId"
                )
                sys.exit(1)
    else:
        if args.referenceId:
            referenceIds = args.referenceId
        else:
            print("No reference id(s) specified with --referenceId.")
            sys.exit(1)

    # We don't have a file of reads, we just want a read filter that we
    # can use to filter the SAM file query sequences.
    reads = parseFASTAFilteringCommandLineOptions(args, Reads())
    samFilter = SAMFilter.parseFilteringOptions(args, reads.filterRead)
    filterAlignment = samFilter.filterAlignment

    proteinGenomeDB = SqliteIndex(args.proteinGenomeDatabase)

    for referenceId in referenceIds:
        try:
            gpInfo = GenomeProteinInfo(
                referenceId, proteinGenomeDB, checkTranslations=args.checkTranslations
            )
        except NoSuchGenomeError:
            print(
                "Reference %r not found in genome database. Ignoring." % referenceId,
                file=sys.stderr,
            )
        else:
            if samfiles:
                if args.progress:
                    print(
                        "Processing %d SAM file%s for matches with %r:"
                        % (
                            len(samfiles),
                            "" if len(samfiles) == 1 else "s",
                            referenceId,
                        ),
                        file=sys.stderr,
                    )
                for i, filename in enumerate(samfiles, start=1):
                    if args.progress:
                        print("  %d: %s" % (i, filename), file=sys.stderr)
                    gpInfo.addSAM(filename, filterAlignment)

            summarize(gpInfo, args.sortOn, args.minReadOffsetCount)

#!/usr/bin/env python

import sys
import warnings
from time import time
from itertools import chain
import argparse

from dark.civ.proteins import SqliteIndexWriter
from dark.taxonomy import (
    addTaxonomyDatabaseCommandLineOptions,
    parseTaxonomyDatabaseCommandLineOptions,
)


def filenamesAndAdders(args, db):
    """
    Get the filenames and database adding functions.

    @param args: A C{Namespace} instance as returned by argparse.
    @param db: An C{SqliteIndexWriter} instance.
    @return: A generator yielding 3-tuples containing a filename, a
        database adder function that can read the file and incorporate it,
        and a C{str} that is either 'gb' or 'json'.
    """
    # Use chain to flatten the lists of lists that we get from using both
    # nargs='+' and action='append'. We use both because it allows people
    # to use (e.g.) --gb on the command line either via "--gb file1 --gb
    # file2" or "--gb file1 file2", or a combination of these. That way
    # it's not necessary to remember which way you're supposed to use it
    # and you also can't be hit by the subtle problem encountered in
    # https://github.com/acorg/dark-matter/issues/453

    if not (args.gb or args.json):
        print(
            "At least one of --gb or --json must be used to indicate input "
            "files to process.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.gb:
        for filename in chain.from_iterable(args.gb):
            yield filename, db.addGenBankFile, "gb"

    if args.json:
        for filename in chain.from_iterable(args.json):
            yield filename, db.addJSONFile, "json"


def main(args, parser):
    """
    Build the protein database.

    @param args: The namespace of command-line arguments returned by
        argparse.parse_args()
    @param parser: An C{argparse.ArgumentParser} instance.
    """

    if (
        args.minGenomeLength is not None
        and args.maxGenomeLength is not None
        and args.minGenomeLength > args.maxGenomeLength
    ):
        raise ValueError("--minGenomeLength cannot be larger than --maxGenomeLength")

    if args.excludeExclusiveHost:
        excludeExclusiveHosts = set(chain.from_iterable(args.excludeExclusiveHost))
    else:
        excludeExclusiveHosts = None

    if args.allowedTaxonomicRanks:
        allowedTaxonomicRanks = set()
        for nameRank in chain.from_iterable(args.allowedTaxonomicRanks):
            name, rank = map(str.lower, map(str.strip, nameRank.split(":", maxsplit=1)))
            allowedTaxonomicRanks.add((name, rank))
    else:
        allowedTaxonomicRanks = None

    taxonomyDatabase = parseTaxonomyDatabaseCommandLineOptions(args, parser)
    progress = args.progress

    if progress:
        overallStart = time()
        totalGenomeCount = totalProteinCount = 0

    with SqliteIndexWriter(args.databaseFile) as db:
        for fileCount, (filename, addFunc, type_) in enumerate(
            filenamesAndAdders(args, db), start=1
        ):
            if args.logFile:
                print("\n>>> Indexing '%s'." % filename, end="\n\n", file=args.logFile)

            if progress:
                start = time()

            examinedGenomeCount, genomeCount, proteinCount = addFunc(
                filename,
                dnaOnly=args.dnaOnly,
                rnaOnly=args.rnaOnly,
                allowedTaxonomicRanks=allowedTaxonomicRanks,
                minGenomeLength=args.minGenomeLength,
                maxGenomeLength=args.maxGenomeLength,
                excludeExclusiveHosts=excludeExclusiveHosts,
                excludeFungusOnlyViruses=args.excludeFungusOnlyViruses,
                excludePlantOnlyViruses=args.excludePlantOnlyViruses,
                databaseName=args.databaseName,
                taxonomyDatabase=taxonomyDatabase,
                proteinSource=args.proteinSource,
                genomeSource=args.genomeSource,
                duplicationPolicy=args.duplicationPolicy,
                logfp=args.logFile,
            )

            if examinedGenomeCount == 0:
                if type_ == "gb":
                    print(
                        "WARNING: No genomes found in %r. Did the GenBank "
                        "download fail on that file?" % filename,
                        file=sys.stderr,
                    )
                else:
                    assert type_ == "json"
                    print(
                        "WARNING: no genomes found in JSON file %r." % filename,
                        file=sys.stderr,
                    )

            if progress:
                elapsed = time() - start
                totalGenomeCount += genomeCount
                totalProteinCount += proteinCount
                print(
                    "Processed %r: added %3d of %3d genome%s (%5d "
                    "protein%s) in %.2f seconds."
                    % (
                        filename,
                        genomeCount,
                        examinedGenomeCount,
                        " " if examinedGenomeCount == 1 else "s",
                        proteinCount,
                        "" if proteinCount == 1 else "s",
                        elapsed,
                    ),
                    file=sys.stderr,
                )

    if progress:
        elapsed = time() - overallStart
        print(
            "%d files (containing %d genomes and %d proteins) "
            "indexed in %.2f seconds (%.2f mins)."
            % (fileCount, totalGenomeCount, totalProteinCount, elapsed, elapsed / 60),
            file=sys.stderr,
        )


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Create a genome/protein sqlite3 database from GenBank "
        "and JSON files. The protein sequences for the database are "
        "printed to standard output for indexing by other tools."
    ),
)

parser.add_argument(
    "--databaseFile",
    required=True,
    help=("The output file. This file must not exist (use --force to " "overwrite)."),
)

parser.add_argument(
    "--databaseName",
    help=(
        "The database that the records in the (--gb) GenBank files came "
        'from (e.g., "refseq" or "RVDB").'
    ),
)

parser.add_argument(
    "--duplicationPolicy",
    choices=("error", "ignore"),
    default="ignore",
    help=(
        "What to do if a to-be-inserted accession number is already "
        'present in the database. "error" results in a ValueError being '
        'raised, "ignore" will ignore the duplicate. It should also be '
        "possible to update (i.e., replace) but that is not supported yet."
    ),
)

parser.add_argument("--progress", action="store_true", help="Print indexing progress.")

parser.add_argument(
    "--logFile",
    type=argparse.FileType("w"),
    help="Write indexing details to a log file.",
)

parser.add_argument(
    "--noWarnings",
    action="store_true",
    help="Do not print warnings about unparseable GenBank or JSON records.",
)

# A mutually exclusive group for DNA/RNA only.
group = parser.add_mutually_exclusive_group()

group.add_argument(
    "--dnaOnly", action="store_true", help="If given, only include DNA viruses."
)

group.add_argument(
    "--rnaOnly", action="store_true", help="If given, only include RNA viruses."
)

parser.add_argument(
    "--maxGenomeLength",
    type=int,
    help="Genomes longer than this will not be considered.",
)

parser.add_argument(
    "--minGenomeLength",
    type=int,
    help="Genomes shorter than this will not be considered.",
)

parser.add_argument(
    "--excludeFungusOnlyViruses",
    action="store_true",
    help=(
        "If given, exclude fungus-only viruses (i.e., viruses that only "
        "infect fungi host species)."
    ),
)

parser.add_argument(
    "--excludePlantOnlyViruses",
    action="store_true",
    help=(
        "If given, exclude plant-only viruses (i.e., viruses that only "
        "infect plant host species)."
    ),
)

parser.add_argument(
    "--gb",
    metavar="GenBank-file",
    nargs="+",
    action="append",
    help=(
        "The GenBank file(s) to make the database from. These may be "
        "uncompressed, or compressed with bgzip (from samtools), with "
        "a .gz suffix."
    ),
)

parser.add_argument(
    "--json",
    metavar="JSON-file",
    nargs="+",
    action="append",
    help=(
        "The JSON file(s) to make the database from. These contain "
        "genome and protein information for cases where sequences that "
        "are not in GenBank should be added."
    ),
)

parser.add_argument(
    "--allowedTaxonomicRanks",
    metavar="name:rank",
    nargs="+",
    action="append",
    help=(
        "Strings of the form name:rank to restrict included viruses to be "
        "from any of a specific set of taxonomic ranks. May be repeated. E.g., "
        "--allowedTaxonomicRanks nidovirales:order "
        "--allowedTaxonomicRanks retroviridae:family "
        "Viruses will be included only if they match at least one of the name, "
        "rank pairs. The strings are case insensitive."
    ),
)

parser.add_argument(
    "--proteinSource",
    default="GenBank",
    help=(
        "The source of the accession numbers for the proteins found in the "
        "input files. This becomes part of the sequence id printed in the "
        "protein FASTA output."
    ),
)

parser.add_argument(
    "--genomeSource",
    default="GenBank",
    help=(
        "The source of the accession numbers for the genomes in the input "
        "files. This becomes part of the sequence id printed in the "
        "protein FASTA output."
    ),
)

parser.add_argument(
    "--excludeExclusiveHost",
    metavar="hostname",
    nargs="+",
    action="append",
    choices=(
        "algae",
        "archaea",
        "bacteria",
        "diatom",
        "environment",
        "eukaryotes",
        "fungi",
        "human",
        "invertebrates",
        "plants",
        "protozoa",
        "vertebrates",
    ),
    help=(
        "A host type that should be excluded, but only if this is the only "
        "host of the virus."
    ),
)

addTaxonomyDatabaseCommandLineOptions(parser)

args = parser.parse_args()

if args.noWarnings:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main(args, parser)
else:
    main(args, parser)

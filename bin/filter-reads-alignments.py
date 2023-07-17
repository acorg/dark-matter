#!/usr/bin/env python

import sys
import argparse
from itertools import chain

from dark.fasta import FastaReads
from dark.fastq import FastqReads


if __name__ == "__main__":
    # We do not use the addFASTACommandLineOptions and
    # parseFASTACommandLineOptions utility functions below because we allow
    # multiple FASTA or FASTQ files on the command line, which we specify
    # by --fasta and --fastq. And those names clash with the option names
    # used by those utility functions.

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Filter reads alignments to select a subset of reads. ",
            "Given BLAST or DIAMOND JSON output files and the "
            "FASTA (or FASTQ) query sequence files given to BLAST or "
            "DIAMOND, and some filtering criteria, print the input "
            "reads that pass the filtering.",
        ),
    )

    parser.add_argument(
        "--matcher",
        default="blast",
        choices=("blast", "diamond"),
        help="The matching algorithm that was used to produce the JSON.",
    )

    parser.add_argument(
        "--json",
        metavar="FILE.JSON[.BZ2]",
        nargs="+",
        action="append",
        required=True,
        help="the JSON file(s) of BLAST or DIAMOND output.",
    )

    # A mutually exclusive group for either FASTA or FASTQ files.
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "--fasta",
        metavar="FILE.FASTA",
        nargs="+",
        action="append",
        help=("the FASTA file(s) of sequences that were given to BLAST " "or DIAMOND."),
    )

    group.add_argument(
        "--fastq",
        metavar="FILE.FASTQ",
        nargs="+",
        action="append",
        help=("the FASTQ file(s) of sequences that were given to BLAST " "or DIAMOND."),
    )

    # Args specific to DIAMOND.

    # A group for either the DIAMOND FASTA file or a sqlite3 database
    # of the FASTA used to make the DIAMOND database.
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--diamondDatabaseFastaFilename",
        metavar="FILE.FASTA",
        help=(
            "The filename of the FASTA file used to make the DIAMOND "
            "database. If --matcher diamond is used, either this argument "
            "or --diamondSqliteDatabaseFilename must be specified."
        ),
    )

    group.add_argument(
        "--diamondSqliteDatabaseFilename",
        metavar="FILE.SQL",
        help=(
            "The filename of the sqlite3 database file of FASTA metadata, "
            "made from the FASTA that was used to make the DIAMOND "
            "database. If --matcher diamond is used, either this argument "
            "or --diamondDatabaseFilename must be specified."
        ),
    )

    parser.add_argument(
        "--diamondDatabaseFastaDirectory",
        metavar="DIR",
        help=(
            "The directory where the FASTA file used to make the DIAMOND "
            "database can be found. This argument is only useful when "
            "--diamondSqliteDatabaseFilename is specified."
        ),
    )

    # Args for filtering on ReadsAlignments.
    parser.add_argument(
        "--minStart",
        type=int,
        metavar="OFFSET",
        help="Reads that start before this subject offset should not be " "shown.",
    )

    parser.add_argument(
        "--maxStop",
        type=int,
        metavar="OFFSET",
        help="Reads that end after this subject offset should not be shown.",
    )

    parser.add_argument(
        "--oneAlignmentPerRead",
        default=False,
        action="store_true",
        help="If True, only keep the best alignment for each read.",
    )

    parser.add_argument(
        "--maxAlignmentsPerRead",
        type=int,
        metavar="N",
        help=(
            "Reads with more than this many alignments will be elided. Pass "
            "zero to only keep reads with no matches (alignments)."
        ),
    )

    parser.add_argument(
        "--scoreCutoff",
        type=float,
        metavar="SCORE",
        help=("A float score. Matches with scores worse than this will be " "ignored."),
    )

    parser.add_argument(
        "--maxHspsPerHit",
        type=int,
        metavar="N",
        help=("A numeric maximum number of HSPs to keep for each subject " "match."),
    )

    parser.add_argument(
        "--whitelist",
        nargs="+",
        action="append",
        metavar="TITLE",
        help="subject titles that should be whitelisted",
    )

    parser.add_argument(
        "--blacklist",
        nargs="+",
        action="append",
        metavar="TITLE",
        help="subject titles that should be blacklisted",
    )

    parser.add_argument(
        "--titleRegex",
        metavar="TITLE-REGEX",
        help="a regex that subject titles must match.",
    )

    parser.add_argument(
        "--negativeTitleRegex",
        metavar="TITLE-REGEX",
        help="a regex that subject titles must not match.",
    )

    parser.add_argument(
        "--truncateTitlesAfter",
        help=(
            "a string that subject titles will be truncated beyond. If the "
            "truncated version of a title has already been seen, "
            "that title will be skipped."
        ),
    )

    parser.add_argument(
        "--minSequenceLen",
        type=int,
        metavar="N",
        help="subjects of lesser length will be elided.",
    )

    parser.add_argument(
        "--maxSequenceLen",
        type=int,
        metavar="N",
        help="subjects of greater length will be elided.",
    )

    parser.add_argument(
        "--taxonomy",
        metavar="NAME",
        help=(
            "the taxonomic group that subjects must match "
            'E.g., "Vira" will filter on viruses.'
        ),
    )

    args = parser.parse_args()

    # Flatten lists of lists that we get from using both nargs='+' and
    # action='append'. We use both because it allows people to use (e.g.)
    # --json on the command line either via "--json file1 --json file2" or
    # "--json file1 file2", or a combination of these. That way it's not
    # necessary to remember which way you're supposed to use it and you also
    # can't be hit by the subtle problem encountered in
    # https://github.com/acorg/dark-matter/issues/453
    jsonFiles = list(chain.from_iterable(args.json))
    whitelist = set(chain.from_iterable(args.whitelist)) if args.whitelist else None
    blacklist = set(chain.from_iterable(args.blacklist)) if args.blacklist else None

    # TODO: Add a --readClass command-line option in case we want to
    # process FASTA containing AA sequences.
    if args.fasta:
        reads = FastaReads(list(chain.from_iterable(args.fasta)))
    else:
        reads = FastqReads(list(chain.from_iterable(args.fastq)))

    if args.matcher == "blast":
        from dark.blast.alignments import BlastReadsAlignments

        readsAlignments = BlastReadsAlignments(reads, jsonFiles)
    else:
        # Must be 'diamond' (due to parser.add_argument 'choices' argument).
        if (
            args.diamondDatabaseFastaFilename is None
            and args.diamondSqliteDatabaseFilename is None
        ):
            print(
                "Either --diamondDatabaseFastaFilename or "
                "--diamondSqliteDatabaseFilename must be used with "
                "--matcher diamond.",
                file=sys.stderr,
            )
            sys.exit(1)
        elif not (
            args.diamondDatabaseFastaFilename is None
            or args.diamondSqliteDatabaseFilename is None
        ):
            print(
                "--diamondDatabaseFastaFilename and "
                "--diamondSqliteDatabaseFilename cannot both be used with "
                "--matcher diamond.",
                file=sys.stderr,
            )
            sys.exit(1)

        from dark.diamond.alignments import DiamondReadsAlignments

        readsAlignments = DiamondReadsAlignments(
            reads,
            jsonFiles,
            databaseFilename=args.diamondDatabaseFastaFilename,
            databaseDirectory=args.diamondDatabaseFastaDirectory,
            sqliteDatabaseFilename=args.diamondSqliteDatabaseFilename,
        )

    readsAlignments.filter(
        maxAlignmentsPerRead=args.maxAlignmentsPerRead,
        minSequenceLen=args.minSequenceLen,
        maxSequenceLen=args.maxSequenceLen,
        minStart=args.minStart,
        maxStop=args.maxStop,
        oneAlignmentPerRead=args.oneAlignmentPerRead,
        maxHspsPerHit=args.maxHspsPerHit,
        scoreCutoff=args.scoreCutoff,
        whitelist=whitelist,
        blacklist=blacklist,
        titleRegex=args.titleRegex,
        negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter,
        taxonomy=args.taxonomy,
    )

    format_ = "fasta" if args.fasta else "fastq"
    write = sys.stdout.write

    for readAlignments in readsAlignments:
        write(readAlignments.read.toString(format_=format_))

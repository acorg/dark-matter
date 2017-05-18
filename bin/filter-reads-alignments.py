#!/usr/bin/env python

"""
Given a BLAST or DIAMOND JSON output files, the corresponding FASTA (or FASTQ)
sequence files, and filtering criteria, filter the reads according to given
criteria and write them to standard output.

Run with --help for help.
"""

from __future__ import print_function

import sys
import argparse
from itertools import chain

from dark.fasta import FastaReads
from dark.fastq import FastqReads


if __name__ == '__main__':

    # We do not use the addFASTACommandLineOptions and
    # parseFASTACommandLineOptions utility functions below because we allow
    # multiple FASTA or FASTQ files on the command line, which we specify
    # by --fasta and --fastq. And those names clash with the option names
    # used by those utility functions.

    parser = argparse.ArgumentParser(
        description='Filter reads alignments to select a subset of reads',
        epilog=('Given BLAST or DIAMOND JSON output files, the '
                'corresponding FASTA (or FASTQ) sequence files, and '
                'filtering  criteria, print a subset of matching reads.'))

    parser.add_argument(
        '--matcher', default='blast', choices=('blast', 'diamond'),
        help='The matching algorithm that was used to produce the JSON.')

    parser.add_argument(
        '--json', metavar='JSON-file', nargs='+', action='append',
        required=True, help='the JSON file(s) of BLAST or DIAMOND output.')

    # A mutually exclusive group for either FASTA or FASTQ files.
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        '--fasta', metavar='FASTA-file', nargs='+', action='append',
        help=('the FASTA file(s) of sequences that were given to BLAST '
              'or DIAMOND.'))

    group.add_argument(
        '--fastq', metavar='FASTQ-file', nargs='+', action='append',
        help=('the FASTQ file(s) of sequences that were given to BLAST '
              'or DIAMOND.'))

    # Args specific to DIAMOND.

    # A group for either the DIAMOND FASTA file or a sqlite3 database
    # of the FASTA used to make the DIAMOND database.
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '--diamondDatabaseFastaFilename',
        help=('The filename of the FASTA file used to make the DIAMOND '
              'database. If --matcher diamond is used, either this argument '
              'or --diamondSqliteDatabaseFilename must be specified.'))

    group.add_argument(
        '--diamondSqliteDatabaseFilename',
        help=('The filename of the sqlite3 database file of FASTA metadata, '
              'made from the FASTA that was used to make the DIAMOND '
              'database. If --matcher diamond is used, either this argument '
              'or --diamondDatabaseFilename must be specified.'))

    parser.add_argument(
        '--diamondDatabaseFastaDirectory',
        help=('The directory where the FASTA file used to make the DIAMOND '
              'database can be found. This argument is only useful when '
              '--diamondSqliteDatabaseFilename is specified.'))

    # Args for filtering on ReadsAlignments.
    parser.add_argument(
        '--minStart', type=int, default=None,
        help='Reads that start before this subject offset should not be '
        'shown.')

    parser.add_argument(
        '--maxStop', type=int, default=None,
        help='Reads that end after this subject offset should not be shown.')

    parser.add_argument(
        '--oneAlignmentPerRead', default=False, action='store_true',
        help='If True, only keep the best alignment for each read.')

    parser.add_argument(
        '--maxAlignmentsPerRead', type=int, default=None,
        help=('Reads with more than this many alignments will be elided. Pass '
              'zero to only keep reads with no matches (alignments).'))

    parser.add_argument(
        '--scoreCutoff', type=float, default=None,
        help=('A float score. Matches with scores worse than this will be '
              'ignored.'))

    parser.add_argument(
        '--maxHspsPerHit', type=int, default=None,
        help='A numeric max number of HSPs to show for each hit on hitId.')

    parser.add_argument(
        '--whitelist', nargs='+', default=None, action='append',
        help='sequence titles that should be whitelisted')

    parser.add_argument(
        '--blacklist', nargs='+', default=None, action='append',
        help='sequence titles that should be blacklisted')

    parser.add_argument(
        '--titleRegex', default=None,
        help='a regex that sequence titles must match.')

    parser.add_argument(
        '--negativeTitleRegex', default=None,
        help='a regex that sequence titles must not match.')

    parser.add_argument(
        '--truncateTitlesAfter', default=None,
        help=('a string that titles will be truncated beyond. If the '
              'truncated version of a title has already been seen, '
              'that title will be skipped.'))

    parser.add_argument(
        '--minSequenceLen', type=int, default=None,
        help='sequences of lesser length will be elided.')

    parser.add_argument(
        '--maxSequenceLen', type=int, default=None,
        help='sequences of greater length will be elided.')

    parser.add_argument(
        '--taxonomy', default=None,
        help=('a string of the taxonomic group on which should be '
              'filtered. eg "Vira" will filter on viruses.'))

    args = parser.parse_args()

    # Flatten lists of lists that we get from using both nargs='+' and
    # action='append'. We use both because it allows people to use (e.g.)
    # --json on the command line either via "--json file1 --json file2" or
    # "--json file1 file2", or a combination of these. That way it's not
    # necessary to remember which way you're supposed to use it and you also
    # can't be hit by the subtle problem encountered in
    # https://github.com/acorg/dark-matter/issues/453
    jsonFiles = list(chain.from_iterable(args.json))
    whitelist = (
        set(chain.from_iterable(args.whitelist)) if args.whitelist else None)
    blacklist = (
        set(chain.from_iterable(args.blacklist)) if args.blacklist else None)

    # TODO: Add a --readClass command-line option in case we want to
    # process FASTA containing AA sequences.
    if args.fasta:
        reads = FastaReads(list(chain.from_iterable(args.fasta)))
    else:
        reads = FastqReads(list(chain.from_iterable(args.fastq)))

    if args.matcher == 'blast':
        from dark.blast.alignments import BlastReadsAlignments
        readsAlignments = BlastReadsAlignments(reads, jsonFiles)
    else:
        # Must be 'diamond' (due to parser.add_argument 'choices' argument).
        if (args.diamondDatabaseFastaFilename is None and
                args.diamondSqliteDatabaseFilename is None):
            print('Either --diamondDatabaseFastaFilename or '
                  '--diamondSqliteDatabaseFilename must be used with '
                  '--matcher diamond.', file=sys.stderr)
            sys.exit(1)
        elif not (args.diamondDatabaseFastaFilename is None or
                  args.diamondSqliteDatabaseFilename is None):
            print('--diamondDatabaseFastaFilename and '
                  '--diamondSqliteDatabaseFilename cannot both be used with '
                  '--matcher diamond.', file=sys.stderr)
            sys.exit(1)

        from dark.diamond.alignments import DiamondReadsAlignments
        readsAlignments = DiamondReadsAlignments(
            reads, jsonFiles,
            databaseFilename=args.diamondDatabaseFastaFilename,
            databaseDirectory=args.diamondDatabaseFastaDirectory,
            sqliteDatabaseFilename=args.diamondSqliteDatabaseFilename)

    readsAlignments.filter(
        maxAlignmentsPerRead=args.maxAlignmentsPerRead,
        minSequenceLen=args.minSequenceLen,
        maxSequenceLen=args.maxSequenceLen,
        minStart=args.minStart, maxStop=args.maxStop,
        oneAlignmentPerRead=args.oneAlignmentPerRead,
        maxHspsPerHit=args.maxHspsPerHit,
        scoreCutoff=args.scoreCutoff,
        whitelist=whitelist, blacklist=blacklist,
        titleRegex=args.titleRegex, negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter, taxonomy=args.taxonomy)

    format_ = 'fasta' if args.fasta else 'fastq'

    for readAlignments in readsAlignments:
        print(readAlignments.read.toStr(format_=format_))

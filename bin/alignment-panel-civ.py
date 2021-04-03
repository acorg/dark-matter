#!/usr/bin/env python

"""
Given a BLAST or DIAMOND JSON output files, the corresponding FASTA (or FASTQ)
sequence files, and filtering criteria, produce a summary of matched titles
and (optionally) an alignment panel.

Run with --help for help.
"""

from __future__ import print_function

import os
import sys
import argparse
from json import dump
from collections import defaultdict
from itertools import chain

# It's not clear that the PDF backend is the right choice here, but it
# works (i.e., the generation of PNG images works fine).
import matplotlib
matplotlib.use('PDF')

# These imports are here because dark.graphics imports matplotlib.pyplot
# and we need to set the matplotlib backend (see above) before that import
# happens. So please don't move these imports higher in this file.
from dark.titles import TitlesAlignments
from dark.fasta import FastaReads
from dark.fastq import FastqReads
from dark.civ.graphics import (
    DEFAULT_LOG_LINEAR_X_AXIS_BASE, alignmentPanelHTML)
from dark.civ.proteins import SqliteIndex
from dark.utils import numericallySortFilenames


def parseColors(colors, args):
    """
    Parse read id color specification.

    @param colors: A C{list} of C{str}s. Each item is of the form, e.g.,
        'green X Y Z...', where each of X, Y, Z, ... etc. is either a read
        id or the name of a FASTA or FASTQ file containing reads whose ids
        should be displayed with the corresponding color. Note that if read
        ids contain spaces you will need to use the latter (i.e. FASTA/Q file
        name) approach because C{args.colors} is split on whitespace.
    @param args: The argparse C{Namespace} instance holding the other parsed
        command line arguments.
    @return: A C{dict} whose keys are colors and whose values are sets of
        read ids.
    """
    result = defaultdict(set)
    for colorInfo in colors:
        readIds = colorInfo.split()
        color = readIds.pop(0)
        for readId in readIds:
            if os.path.isfile(readId):
                filename = readId
                if args.fasta:
                    reads = FastaReads(filename)
                else:
                    reads = FastqReads(filename)
                for read in reads:
                    result[color].add(read.id)
            else:
                result[color].add(readId)
    return result


if __name__ == '__main__':

    # We do not use the addFASTACommandLineOptions and
    # parseFASTACommandLineOptions utility functions below because we allow
    # multiple FASTA or FASTQ files on the command line, which we specify
    # by --fasta and --fastq. And those names clash with the option names
    # used by those utility functions.

    parser = argparse.ArgumentParser(
        description='Produce an alignment panel in HTML',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=('Given BLAST or DIAMOND JSON output files, the '
                'corresponding FASTA (or FASTQ) sequence files, and '
                'filtering criteria, produce an alignment panel.'))

    parser.add_argument(
        '--proteinGenomeDatabase', required=True,
        help=('The filename of an Sqlite3 database holding protein and '
              'genome information, as built by make-protein-database.py'))

    parser.add_argument(
        '--earlyExit', default=False, action='store_true',
        help=('If True, just print the number of interesting matches, but do '
              'not create the alignment panel.'))

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
        help=('the FASTA file(s) of read sequences that were given to BLAST '
              'or DIAMOND.'))

    group.add_argument(
        '--fastq', metavar='FASTQ-file', nargs='+', action='append',
        help=('the FASTQ file(s) of read sequences that were given to BLAST '
              'or DIAMOND.'))

    # Args for filtering on ReadsAlignments.
    parser.add_argument(
        '--minStart', type=int,
        help='Reads that start before this subject offset should not be '
        'shown.')

    parser.add_argument(
        '--maxStop', type=int,
        help='Reads that end after this subject offset should not be shown.')

    parser.add_argument(
        '--oneAlignmentPerRead', default=False, action='store_true',
        help='If True, only keep the best alignment for each read.')

    parser.add_argument(
        '--maxAlignmentsPerRead', type=int,
        help=('Reads with more than this many alignments will be elided. Pass '
              'zero to only keep reads with no matches (alignments).'))

    parser.add_argument(
        '--scoreCutoff', type=float,
        help=('A float score. Matches with scores worse than this will be '
              'ignored.'))

    parser.add_argument(
        '--percentageIdenticalCutoff', type=float,
        help=('A float percentage identity (0.0 to 100.0). Matches with '
              'percent identity scores less than this will be ignored.'))

    parser.add_argument(
        '--percentagePositiveCutoff', type=float,
        help=('A float percentage (0.0 to 100.0). Matches with percent '
              'positive scores less than this will be ignored. A single AA '
              'match is considered positive (by DIAMOND) if its entry in the '
              'substitution (scoring) matrix is positive. DIAMOND uses '
              'the BLOSUM62 matrix by default.'))

    parser.add_argument(
        '--maxHspsPerHit', type=int,
        help='A numeric max number of HSPs to show for each hit on hitId.')

    parser.add_argument(
        '--whitelist', nargs='+', action='append',
        help='sequence titles that should be whitelisted')

    parser.add_argument(
        '--whitelistFile',
        help=('A file containing sequence titles (one per line) to be '
              'whitelisted'))

    parser.add_argument(
        '--blacklist', nargs='+', action='append',
        help='sequence titles that should be blacklisted')

    parser.add_argument(
        '--blacklistFile',
        help=('A file containing sequence titles (one per line) to be '
              'blacklisted'))

    parser.add_argument(
        '--titleRegex',
        help='a regex that sequence titles must match.')

    parser.add_argument(
        '--negativeTitleRegex',
        help='a regex that sequence titles must not match.')

    parser.add_argument(
        '--truncateTitlesAfter',
        help=('a string that titles will be truncated beyond. If the '
              'truncated version of a title has already been seen, '
              'that title will be skipped.'))

    parser.add_argument(
        '--minSequenceLen', type=int,
        help='sequences of lesser length will be elided.')

    parser.add_argument(
        '--maxSequenceLen', type=int,
        help='sequences of greater length will be elided.')

    parser.add_argument(
        '--taxonomy',
        help=('a string of the taxonomic group on which should be '
              'filtered. eg "Vira" will filter on viruses.'))

    # Args for filtering on TitlesAlignments.
    parser.add_argument(
        '--minMatchingReads', type=int,
        help='sequences that are matched by fewer reads will be elided.')

    parser.add_argument(
        '--maxMatchingReads', type=int,
        help='sequences that are matched by more reads will be elided.')

    parser.add_argument(
        '--minMedianScore', type=float,
        help=('sequences that are matched with a median score that is '
              'worse will be elided.'))

    parser.add_argument(
        '--withScoreBetterThan', type=float,
        help=('sequences that are matched without at least one score '
              'at least this good will be elided.'))

    parser.add_argument(
        '--minNewReads', type=float,
        help=('The fraction of its reads by which a new read set must differ '
              'from all previously seen read sets in order to be considered '
              'acceptably different.'))

    parser.add_argument(
        '--maxTitles', type=int,
        help=('The maximum number of titles to keep. If more titles than '
              'this result from the filtering, titles will be sorted '
              '(according to the --sortOn value) and only the best will be '
              'retained.'))

    parser.add_argument(
        '--minCoverage', type=float,
        help=('The (0.0 to 1.0) minimum fraction of a subject sequence that '
              'must be matched by at least one read.'))

    # Args for the alignment panel
    parser.add_argument(
        '--sortOn', default='maxScore',
        choices=('maxScore', 'medianScore', 'readCount', 'length', 'title'),
        help='The attribute to sort subplots on.')

    parser.add_argument(
        '--rankValues', type=bool, default=False,
        help=('If True, display reads with a Y axis coord that is the rank of '
              'the score.'))

    parser.add_argument(
        '--outputDir', required=True,
        help='Specifies a directory to write the HTML summary to. Required.')

    parser.add_argument(
        '--color', action='append',
        help=('a string which has a color as the first element and read ids '
              'and/or FASTA file names as the following element(s), separated '
              'by spaces.'))

    parser.add_argument(
        '--equalizeXAxes', default=False, action='store_true',
        help=('If True, all alignment graphs will have their X axes drawn '
              'with the same range.'))

    parser.add_argument(
        '--xRange', default='subject',
        choices=('reads', 'subject'),
        help=('Set the X axis range to show either the subject or the extent '
              'of the reads that hit the subject.'))

    parser.add_argument(
        '--logLinearXAxis', default=False, action='store_true',
        help=('If True, convert read offsets so that empty regions in the '
              'alignment panel plots will only be as wide as their logged '
              'actual values'))

    parser.add_argument(
        '--logBase', type=float, default=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
        help='The base of the logarithm to use if logLinearXAxis is True')

    parser.add_argument(
        '--showFeatures', default=False, action='store_true',
        help=('If specified, look up features for the individual images in '
              'the alignment panel.'))

    parser.add_argument(
        '--sortFilenames', default=False, action='store_true',
        help=('If specified, the JSON and FASTA/Q file names will be '
              'processed in sorted order. The sorting is based on finding '
              'a numerical prefix in the filename. This can be useful when '
              'processing output files produced by systems like HTCondor, '
              'which makes files with names like 1.out, 10.out, etc. that do '
              'not sort properly and so cannot conveniently be given to this '
              'program along with a single FASTA/Q file (because the order of '
              'the results in the files from HTCondor does not match the '
              'order of sequences in the FASTA/Q file.'))

    parser.add_argument(
        '--titlesJSONFile',
        help=('Give a file name for JSON holding information about titles to '
              'be written to'))

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
        if args.sortFilenames:
            files = numericallySortFilenames(chain.from_iterable(args.fasta))
        else:
            files = list(chain.from_iterable(args.fasta))
        reads = FastaReads(files)
    else:
        if args.sortFilenames:
            files = numericallySortFilenames(chain.from_iterable(args.fastq))
        else:
            files = list(chain.from_iterable(args.fastq))
        reads = FastqReads(files)

    if args.matcher == 'blast':
        from dark.blast.alignments import BlastReadsAlignments
        readsAlignments = BlastReadsAlignments(
            reads, jsonFiles,
            sqliteDatabaseFilename=args.sqliteDatabaseFilename,
            sortBlastFilenames=args.sortFilenames)
    else:
        # Must be 'diamond' (due to parser.add_argument 'choices' argument).
        assert args.matcher == 'diamond'
        from dark.diamond.alignments import DiamondReadsAlignments
        readsAlignments = DiamondReadsAlignments(
            reads, jsonFiles, sortFilenames=args.sortFilenames)

    readsAlignments.filter(
        maxAlignmentsPerRead=args.maxAlignmentsPerRead,
        minSequenceLen=args.minSequenceLen,
        maxSequenceLen=args.maxSequenceLen,
        minStart=args.minStart, maxStop=args.maxStop,
        oneAlignmentPerRead=args.oneAlignmentPerRead,
        maxHspsPerHit=args.maxHspsPerHit,
        scoreCutoff=args.scoreCutoff,
        percentageIdenticalCutoff=args.percentageIdenticalCutoff,
        percentagePositiveCutoff=args.percentagePositiveCutoff,
        whitelist=whitelist, blacklist=blacklist,
        whitelistFile=args.whitelistFile, blacklistFile=args.blacklistFile,
        titleRegex=args.titleRegex, negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter, taxonomy=args.taxonomy)

    titlesAlignments = TitlesAlignments(readsAlignments).filter(
        minMatchingReads=args.minMatchingReads,
        maxMatchingReads=args.maxMatchingReads,
        minMedianScore=args.minMedianScore,
        withScoreBetterThan=args.withScoreBetterThan,
        minNewReads=args.minNewReads, maxTitles=args.maxTitles,
        sortOn=args.sortOn, minCoverage=args.minCoverage)

    nTitles = len(titlesAlignments)
    print('Found %d interesting title%s.' %
          (nTitles, '' if nTitles == 1 else 's'), file=sys.stderr)

    # If we've been asked to save the matched title information in JSON,
    # write it out.
    if args.titlesJSONFile:
        with open(args.titlesJSONFile, 'w') as fp:
            dump(titlesAlignments.toDict(), fp, sort_keys=True, indent=2,
                 separators=(',', ': '))

    if nTitles:
        # Use flush=True on the print of the titles so the output is
        # definitely written out. This is because in some cases the
        # following aligment panel generation takes a long time and exceeds
        # its CPU or memory limit on an HPCS system.  Also, it can be good
        # to see the matching title details while waiting for the panel
        # plots to be generated.
        print(titlesAlignments.tabSeparatedSummary(sortOn=args.sortOn),
              flush=True)

    if args.earlyExit:
        sys.exit(0)

    if nTitles == 0:
        print('No alignment panel generated due to no matching titles.',
              file=sys.stderr)
        sys.exit(0)

    idList = parseColors(args.color, args) if args.color else None

    alignmentPanelHTML(
        titlesAlignments, SqliteIndex(args.proteinGenomeDatabase),
        sortOn=args.sortOn, outputDir=args.outputDir,
        idList=idList, equalizeXAxes=args.equalizeXAxes, xRange=args.xRange,
        logLinearXAxis=args.logLinearXAxis, logBase=args.logBase,
        showFeatures=args.showFeatures)

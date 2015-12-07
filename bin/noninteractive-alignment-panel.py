#!/usr/bin/env python

"""
Given a JSON BLAST output file, a FASTA sequence file, and interesting
criteria, produce an alignment panel.

Run with --help for help.
"""

from __future__ import print_function

import os
import sys
import argparse
from collections import defaultdict

from dark.fasta import FastaReads
from dark.blast.alignments import BlastReadsAlignments
from dark.titles import TitlesAlignments
from dark.graphics import DEFAULT_LOG_LINEAR_X_AXIS_BASE, alignmentPanel


def parseColors(colors):
    """
    Parse read id color specification.

    @param colors: A C{list} of C{str}s. Each item is of the form, e.g.,
        'green X Y Z...', where each of X, Y, Z, ... etc. is either a read
        id or the name of a file containing read identifiers one per line.
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
                for read in FastaReads(filename):
                    result[color].add(read.id)
            else:
                result[color].add(readId)
    return result


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Non-interactively generate an alignment panel',
        epilog='Given a FASTA sequence file, a JSON BLAST output files, '
        'and filtering criteria, produce an alignment panel.'
    )

    parser.add_argument(
        '--earlyExit', default=False, action='store_true',
        help='If True, just print the number of interesting matches, but do '
        'not create the alignment panel.')

    # Args for the JSON BLAST and FASTA files.
    parser.add_argument(
        'json', metavar='BLAST-JSON-file', type=str, nargs='+',
        help='the JSON file of BLAST output.')

    parser.add_argument(
        '--fasta', metavar='FASTA-file', type=str, required=True,
        help='the FASTA file of sequences that were given to BLAST.')

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
        help='If C{True}, only keep the best alignment for each read.')

    parser.add_argument(
        '--scoreCutoff', type=float, default=None,
        help=('A float score. Matches with scores worse than '
              'this will be ignored.'))

    parser.add_argument(
        '--maxHspsPerHit', type=int, default=None,
        help='A numeric max number of HSPs to show for each hit on hitId.')

    parser.add_argument(
        '--whitelist', type=str, nargs='+', default=None,
        help='sequence titles that should be whitelisted')

    parser.add_argument(
        '--blacklist', type=str, nargs='+', default=None,
        help='sequence titles that should be blacklisted')

    parser.add_argument(
        '--titleRegex', type=str, default=None,
        help='a regex that sequence titles must match.')

    parser.add_argument(
        '--negativeTitleRegex', type=str, default=None,
        help='a regex that sequence titles must not match.')

    parser.add_argument(
        '--truncateTitlesAfter', type=str, default=None,
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
        '--taxonomy', type=str, default=None,
        help='a string of the taxonomic group on which should be '
        'filtered. eg "Vira" will filter on viruses.')

    # Args for filtering on TitlesAlignments.
    parser.add_argument(
        '--minMatchingReads', type=int, default=None,
        help='sequences that are matched by fewer reads will be elided.')

    parser.add_argument(
        '--minMedianScore', type=float, default=None,
        help='sequences that are matched with a median score that is '
        'worse will be elided.')

    parser.add_argument(
        '--withScoreBetterThan', type=float, default=None,
        help='sequences that are matched without at least one score '
        'at least this good will be elided.')

    parser.add_argument(
        '--minNewReads', type=float, default=None,
        help='The fraction of its reads by which a new read set must differ '
        'from all previously seen read sets in order to be considered '
        'acceptably different.')

    # Args for the alignment panel
    parser.add_argument(
        '--sortOn', type=str, default='maxScore',
        choices=['maxScore', 'medianScore', 'readCount', 'length', 'title'],
        help='The attribute to sort subplots on.')

    parser.add_argument(
        '--rankValues', type=bool, default=False,
        help='If True, display reads with a Y axis coord that is the rank of '
        'the score.')

    parser.add_argument(
        '--outputDir', type=str, default=None,
        help='Specifies a directory to write the HTML summary to.')

    parser.add_argument(
        '--color', type=str, action='append',
        help='a string which has a color as the first element and readIds '
        'or a fastafile as the following element(s), separated by spaces.')

    parser.add_argument(
        '--equalizeXAxes', default=False, action='store_true',
        help='If True, all alignment graphs will have their X axes drawn with '
        'the same range.')

    parser.add_argument(
        '--xRange', type=str, default='subject',
        choices=['reads', 'subject'],
        help='Set the X axis range to show either the subject or the extent '
        'of the reads that hit the subject.')

    parser.add_argument(
        '--logLinearXAxis', default=False, action='store_true',
        help='If True, convert read offsets so that empty regions in the '
        'alignment panel plots will only be as wide as their logged actual '
        'values')

    parser.add_argument(
        '--logBase', type=float, default=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
        help='The base of the logarithm to use if logLinearXAxis is True')

    args = parser.parse_args()
    readsAlignments = BlastReadsAlignments(FastaReads(args.fasta), args.json)

    readsAlignments.filter(
        minSequenceLen=args.minSequenceLen,
        maxSequenceLen=args.maxSequenceLen,
        minStart=args.minStart,
        maxStop=args.maxStop,
        oneAlignmentPerRead=args.oneAlignmentPerRead,
        maxHspsPerHit=args.maxHspsPerHit,
        scoreCutoff=args.scoreCutoff,
        whitelist=set(args.whitelist) if args.whitelist else None,
        blacklist=set(args.blacklist) if args.blacklist else None,
        titleRegex=args.titleRegex,
        negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter,
        taxonomy=args.taxonomy)

    titlesAlignments = TitlesAlignments(readsAlignments).filter(
        minMatchingReads=args.minMatchingReads,
        minMedianScore=args.minMedianScore,
        withScoreBetterThan=args.withScoreBetterThan,
        minNewReads=args.minNewReads)

    nTitles = len(titlesAlignments)
    print('Found %d interesting title%s.' % (nTitles,
                                             '' if nTitles == 1 else 's'))

    if args.earlyExit:
        print('Matched titles (sorted by best score, descending):')
        print('\n'.join(titlesAlignments.sortTitles('maxScore')))
        sys.exit(0)

    alignmentPanel(titlesAlignments, sortOn=args.sortOn, interactive=True,
                   outputDir=args.outputDir,
                   idList=parseColors(args.color) if args.color else None,
                   equalizeXAxes=args.equalizeXAxes, xRange=args.xRange,
                   logLinearXAxis=args.logLinearXAxis, logBase=args.logBase)

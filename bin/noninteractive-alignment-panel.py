#!/usr/bin/env python

"""
Given a JSON BLAST output file, a FASTA sequence file, and interesting
criteria, produce an alignment panel.

Run with --help for help.
"""


def parseStringToIdList(idList):
    """
    Parses the string given in --idList to a dict.

    @param idList: what is specified by the --idList parameter
    """
    idDict = {}
    index = 0
    splittedIdList = idList.split('|')
    if len(splittedIdList) % 2 != 0:
        print >>sys.stderr, 'idList contains %d elements' % len(splittedIdList)
    else:
        while index < len(splittedIdList) - 1:
            color = splittedIdList[index]
            index += 1
            values = splittedIdList[index]
            index += 1
            value = values.split(' ')
            reads = []
            for read in value:
                reads.append(read)
            idDict[color] = reads
    for color, reads in idDict.iteritems():
        if os.path.isfile(reads[0]):
            readIds = []
            with open(reads[0], 'r') as fp:
                for line in fp:
                    readId = line.rstrip()
                    readIds.append(readId)
            idDict[color] = readIds
    return idDict


if __name__ == '__main__':
    import sys
    from dark.graphics import alignmentPanel, report
    from dark.blast import BlastRecords
    import argparse

    parser = argparse.ArgumentParser(
        description='Non-interactively generate an alignment panel',
        epilog='Given a JSON BLAST output file, a FASTA sequence file, '
        'and interesting criteria, produce an alignment panel.'
    )

    # Args for the JSON BLAST and FASTA files.
    parser.add_argument(
        'json', metavar='BLAST-JSON-file', type=str, nargs='+',
        help='the JSON file of BLAST output.')

    parser.add_argument(
        '--fasta', metavar='FASTA-file', type=str,
        help='the FASTA file of sequences that were given to BLAST.')

    # Args for the interesting selection.
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
        '--minSequenceLen', type=int, default=None,
        help='sequences of lesser length will be elided.')

    parser.add_argument(
        '--maxSequenceLen', type=int, default=None,
        help='sequences of greater length will be elided.')

    parser.add_argument(
        '--minMatchingReads', type=int, default=None,
        help='sequences that are matched by fewer reads will be elided.')

    parser.add_argument(
        '--maxMeanEValue', type=float, default=None,
        help='sequences that are matched with a mean e-value that is '
        'greater will be elided.')

    parser.add_argument(
        '--maxMedianEValue', type=float, default=None,
        help='sequences that are matched with a median e-value that is '
        'greater will be elided.')

    parser.add_argument(
        '--withEBetterThan', type=float, default=None,
        help='sequences that are matched without at least one e-value '
        'at least this good will be elided.')

    parser.add_argument(
        '--minMeanBitScore', type=float, default=None,
        help='sequences that are matched with a mean bit score that is '
        'less will be elided.')

    parser.add_argument(
        '--minMedianBitScore', type=float, default=None,
        help='sequences that are matched with a median bit score that is '
        'less will be elided.')

    parser.add_argument(
        '--withBitScoreBetterThan', type=float, default=None,
        help='sequences that are matched without at least one bit score '
        'at least this good will be elided.')

    parser.add_argument(
        '--negativeTitleRegex', type=str, default=None,
        help='a regex that sequence titles must not match.')

    parser.add_argument(
        '--truncateTitlesAfter', type=str, default=None,
        help='a string that titles will be truncated beyond. If a truncated '
        'title has already been seen, that title will be skipped.')

    parser.add_argument(
        '--minNewReads', type=float, default=None,
        help='The fraction of its reads by which a new read set must differ '
        'from all previously seen read sets in order to be considered '
        'acceptably different.')

    # Args for the alignment panel
    parser.add_argument(
        '--db', type=str, default='nt', help='the BLAST db that was used')

    parser.add_argument(
        '--eCutoff', type=float, default=None,
        help='Ignore hits with e-values greater (i.e., worse) than or equal '
        'to this.')

    parser.add_argument(
        '--maxHspsPerHit', type=int, default=None,
        help='A numeric max number of HSPs to show for each hit on hitId.')

    parser.add_argument(
        '--minStart', type=int, default=None,
        help='Reads that start before this subject offset should not be '
        'shown.')

    parser.add_argument(
        '--maxStop', type=int, default=None,
        help='Reads that end after this subject offset should not be shown.')

    parser.add_argument(
        '--sortOn', type=str, default='bitScoreMax',
        choices=['eMean', 'eMedian', 'eMin', 'bitScoreMax', 'bitScoreMean',
                 'bitScoreMedian', 'readCount', 'length', 'title'],
        help='The attribute to sort subplots on.')

    parser.add_argument(
        '--rankValues', type=bool, default=False,
        help='If True, display reads with a Y axis coord that is the rank of '
        'the e-value or bit score.')

    parser.add_argument(
        '--outputDir', type=str, default='.',
        help='Specifies a directory to write the HTML summary to.')

    parser.add_argument(
        '--idList', type=str, default=False,
        help='string of the format "green|readId1 readId2|red|filename.txt". '
        'Allows reading readIds specified directly on the command line, '
        'as well as in a file. If a file is used, each readId must be '
        'followed by a new line.')

    parser.add_argument(
        '--earlyExit', default=False, action='store_true',
        help='If True, just print the number of interesting hits, but do not '
        'create the alignment panel.')

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
        '--plot', type=str, default='e values',
        choices=['bit scores', 'e values'],
        help='What to plot on the Y axis of alignment graphs.')

    args = parser.parse_args()

    blastRecords = BlastRecords(args.json, fastaFilename=args.fasta,
                                blastDb=args.db)

    # Convert white/blacklists lists to sets.
    if args.whitelist is not None:
        args.whitelist = set(args.whitelist)
    if args.blacklist is not None:
        args.blacklist = set(args.blacklist)

    hits = blastRecords.filterHits(
        whitelist=args.whitelist,
        blacklist=args.blacklist,
        minSequenceLen=args.minSequenceLen,
        maxSequenceLen=args.maxSequenceLen,
        minMatchingReads=args.minMatchingReads,
        maxMeanEValue=args.maxMeanEValue,
        maxMedianEValue=args.maxMedianEValue,
        withEBetterThan=args.withEBetterThan,
        minMeanBitScore=args.minMeanBitScore,
        minMedianBitScore=args.minMedianBitScore,
        withBitScoreBetterThan=args.withBitScoreBetterThan,
        titleRegex=args.titleRegex,
        negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter,
        minNewReads=args.minNewReads)

    nHits = len(hits)
    if nHits == 0:
        print >>sys.stderr, 'No interesting hits found. Relax your search!'
        sys.exit(0)

    report('Found %d interesting hit%s.' % (nHits, '' if nHits == 1 else 's'))

    if args.earlyExit:
        print 'Hit titles (sorted by best e-value, descending):'
        print '\n'.join(hits.sortTitles('eMin'))
        sys.exit(0)

    hits.computePlotInfo(
        eCutoff=args.eCutoff, maxHspsPerHit=args.maxHspsPerHit,
        minStart=args.minStart, maxStop=args.maxStop,
        rankValues=args.rankValues)

    if args.idList:
        import os
        args.idList = parseStringToIdList(args.idList)

    alignmentPanel(hits, sortOn=args.sortOn, interactive=True,
                   outputDir=args.outputDir, idList=args.idList,
                   equalizeXAxes=args.equalizeXAxes, xRange=args.xRange,
                   plot=args.plot)

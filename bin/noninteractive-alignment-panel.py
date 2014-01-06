#!/usr/bin/env python

"""
Given a JSON BLAST output file, a FASTA sequence file, and interesting
criteria, produce an alignment panel.

Run with --help for help.
"""

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
        'json', metavar='BLAST-JSON-file', type=str,
        help='the JSON file of BLAST output.')

    parser.add_argument(
        'fasta', metavar='FASTA-file', type=str,
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
        '--maxMinEValue', type=float, default=None,
        help='sequences that are matched with a minimum e-value that is '
        'greater will be elided.')

    parser.add_argument(
        '--negativeTitleRegex', type=str, default=None,
        help='a regex that sequence titles must not match.')

    parser.add_argument(
        '--truncateTitlesAfter', type=str, default=None,
        help='a string that titles will be truncated beyond. If a truncated '
        'title has already been seen, that title will be skipped.')

    # Args for the alignment panel
    parser.add_argument(
        '--db', type=str, default='nt', help='the BLAST db that was used')

    parser.add_argument(
        '--eCutoff', type=float, default=2.0,
        help='converted e values less than this will be ignored.')

    parser.add_argument(
        '--maxHspsPerHit', type=str, default=None,
        help='A numeric max number of HSPs to show for each hit on hitId.')

    parser.add_argument(
        '--minStart', type=str, default=None,
        help='Reads that start before this subject offset should not be '
        'shown.')

    parser.add_argument(
        '--maxStop', type=str, default=None,
        help='Reads that end after this subject offset should not be shown.')

    parser.add_argument(
        '--sortOn', type=str, default='eMedian',
        choices=['eMean', 'eMedia', 'title', 'reads'],
        help='The attribute to sort subplots on.')

    parser.add_argument(
        '--rankEValues', type=bool, default=False,
        help='If True, display reads with a Y axis coord that is the rank of '
        'the e value (sorted decreasingly).')

    parser.add_argument(
        '--outputDir', type=str, default='.',
        help='Specifies a directory to write the HTML summary to.')

    parser.add_argument(
        '--idList', type=str, default=False,
        help='a dictionary. The keys is a color and the values is a list of '
        'read identifiers that should be colored in the respective color.')

    parser.add_argument(
        '--earlyExit', type=bool, default=False, const=True, nargs='?',
        help='If True, just print the number of interesting hits, but do not '
        'create the alignment panel.')

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
        maxMinEValue=args.maxMinEValue,
        titleRegex=args.titleRegex,
        negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter)

    hits.computePlotInfo(
        eCutoff=args.eCutoff, maxHspsPerHit=args.maxHspsPerHit,
        minStart=args.minStart, maxStop=args.maxStop,
        rankEValues=args.rankEValues)

    nInteresting = len(hits)
    if nInteresting == 0:
        print >>sys.stderr, 'No interesting hits found. Relax your search!'
        sys.exit(0)

    report('Found %d interesting hit%s.' %
           (nInteresting, '' if nInteresting == 1 else 's'))

    if args.earlyExit:
        print 'Hit titles:'
        print '\n'.join(sorted(hits.titles.keys()))
        sys.exit(0)

    alignmentPanel(hits, idList=args.idList, interactive=True,
                   outputDir=args.outputDir, sortOn=args.sortOn)

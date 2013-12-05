#!/usr/bin/env python

"""
Given a JSON BLAST output file, a FASTA sequence file, and interesting
criteria, produce an alignment panel.

Run with --help for help.
"""

if __name__ == '__main__':
    import sys
    from Bio import SeqIO
    import os.path
    from json import loads, dumps
    from dark.utils import (alignmentPanel, getAllHitsForSummary,
                            interestingRecords, report, summarizeAllRecords)
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

    # Args for summarizeReads.
    parser.add_argument(
        '--eCutoff_for_summary', type=float, default=None,
        help='eCutoff passed to summarizeAllRecords.')

    # Args for the interesting selection.
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
        '--negativeTitleRegex', type=str, default=None,
        help='a regex that sequence titles must not match.')

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

    parser.add_argument(
        '--summaryFile', type=str, default=None,
        help='A file to write the summary from summarizeAllRecords to.')

    args = parser.parse_args()
    report('Reading FASTA from %r.' % args.fasta)
    fasta = list(SeqIO.parse(args.fasta, 'fasta'))
    report('Read %d sequences.' % len(fasta))
    report('Summarizing all records.')

    if args.summaryFile:
        if os.path.isfile(args.summaryFile):
            with open(args.summaryFile, 'r') as f:
                summaryString = f.read()
                summary = loads(summaryString)
        else:
            summary = summarizeAllRecords(args.json, args.eCutoff_for_summary)
            with open(args.summaryFile, 'w') as f:
                f.write(dumps(summary))
    else:
        summary = summarizeAllRecords(args.json, args.eCutoff_for_summary)

    interesting = interestingRecords(
        summary, titleRegex=args.titleRegex,
        minSequenceLen=args.minSequenceLen, maxSequenceLen=args.maxSequenceLen,
        minMatchingReads=args.minMatchingReads,
        maxMeanEValue=args.maxMeanEValue, maxMedianEValue=args.maxMedianEValue,
        negativeTitleRegex=args.negativeTitleRegex)

    nInteresting = len(interesting)
    if nInteresting == 0:
        print >>sys.stderr, 'No interesting hits found. Relax your search!'
        sys.exit(0)

    report('Found %d interesting hit%s.' %
           (nInteresting, '' if nInteresting == 1 else 's'))

    if args.earlyExit:
        sys.exit(0)

    allHits = getAllHitsForSummary(interesting, args.json)
    alignmentPanel(
        interesting, allHits, fasta, db=args.db, eCutoff=args.eCutoff,
        maxHspsPerHit=args.maxHspsPerHit, minStart=args.minStart,
        maxStop=args.maxStop, sortOn=args.sortOn, rankEValues=args.rankEValues,
        interactive=True, outputDir=args.outputDir, idList=args.idList)

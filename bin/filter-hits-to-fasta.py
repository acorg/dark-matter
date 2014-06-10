#!/usr/bin/env python

"""
Given a JSON BLAST output file, a FASTA sequence file, and interesting
criteria, write a FASTA file (to sys.stdout) of the reads that match the
criteria.

Run with --help for help.
"""

if __name__ == '__main__':
    from Bio import SeqIO
    from dark.blast import BlastRecords
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='Extract FASTA from matching BLAST hits.',
        epilog='Given a JSON BLAST output file, a FASTA sequence file, '
        'and interesting criteria, produce a FASTA file on stdout of '
        'the reads that match the criteria.'
    )

    # Args for the JSON BLAST and FASTA files.
    parser.add_argument(
        'json', metavar='BLAST-JSON-file', type=str, nargs='+',
        help='the JSON file of BLAST output.')

    parser.add_argument(
        '--fasta', metavar='FASTA-file', type=str, required=True,
        help='the FASTA file of sequences that were given to BLAST.')

    parser.add_argument(
        '--db', type=str, default='nt', help='the BLAST db that was used')

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

    # Args for computePlotInfo
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

    readNums = set()

    if len(hits):
        # If we've been given any relevant arguments, call computePlotInfo
        # to do further filtering. Then pull all read numbers from either
        # the hitInfo or the plotInfo (depending on whether computePlotInfo
        # was called).  This is a speed-up to avoid calling computePlotInfo
        # unless we need to (as it re-reads all records).
        if (args.eCutoff is None and args.maxHspsPerHit is None and
                args.minStart is None and args.maxStop is None):
            # No need to call computePlotInfo, all read numbers are already in
            # the hitInfo for each title.
            for hitInfo in hits.titles.itervalues():
                readNums.update(hitInfo['readNums'])
        else:
            hits.computePlotInfo(
                eCutoff=args.eCutoff, maxHspsPerHit=args.maxHspsPerHit,
                minStart=args.minStart, maxStop=args.maxStop)
            for hitInfo in hits.titles.itervalues():
                plotInfo = hitInfo['plotInfo']
                if plotInfo is not None:
                    readNums.update(plotInfo['readNums'])

    # Write a FASTA file for all the matched read numbers. Note that read
    # numbers start at zero.
    found = []
    with open(args.fasta) as fp:
        for readNum, seq in enumerate(SeqIO.parse(fp, 'fasta')):
            if readNum in readNums:
                readNums.remove(readNum)
                found.append(seq)
                # Break out of the loop early if we've already found all
                # required reads.
                if not readNums:
                    break

    # If any read numbers are left in readNums, we have a problem.
    if readNums:
        print >>sys.stderr, (
            'WARNING: %d read number%s not found. It is very likely you have '
            'given a FASTA file (%s) that was not the one used to generate '
            'the BLAST output in %s. Read numbers not found: %s' %
            (len(readNums),
             '' if len(readNums) == 1 else 's were',
             args.fasta,
             args.json,
             ' '.join(sorted(readNums))))
        sys.exit(1)
    else:
        SeqIO.write(found, sys.stdout, 'fasta')
        print >>sys.stderr, 'Found %d matching reads.' % len(found)

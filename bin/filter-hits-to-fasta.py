#!/usr/bin/env python

"""
Given a JSON BLAST output file, a FASTA sequence file, and interesting
criteria, write a FASTA file (to sys.stdout) of the reads that match the
criteria.

Run with --help for help.
"""

from __future__ import print_function

import sys
import argparse

from dark.reads import Reads
from dark.fasta import FastaReads
from dark.titles import TitlesAlignments
from dark.blast.alignments import BlastReadsAlignments


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Extract FASTA from matching BLAST hits.',
        epilog='Given JSON BLAST output files, a FASTA sequence file, '
        'and filtering criteria, produce a FASTA file on stdout '
        'containing the reads that match the criteria.'
    )

    # Args for the JSON BLAST and FASTA files.
    parser.add_argument(
        'json', metavar='BLAST-JSON-file', nargs='+',
        help='the JSON file of BLAST output.')

    parser.add_argument(
        '--fasta', metavar='FASTA-file', required=True,
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
        '--whitelist', nargs='+', default=None,
        help='sequence titles that should be whitelisted')

    parser.add_argument(
        '--blacklist', nargs='+', default=None,
        help='sequence titles that should be blacklisted')

    parser.add_argument(
        '--titleRegex', default=None,
        help='a regex that sequence titles must match.')

    parser.add_argument(
        '--negativeTitleRegex', default=None,
        help='a regex that sequence titles must not match.')

    parser.add_argument(
        '--truncateTitlesAfter', default=None,
        help='a string that titles will be truncated beyond. If a truncated '
        'title has already been seen, that title will be skipped.')

    parser.add_argument(
        '--minSequenceLen', type=int, default=None,
        help='sequences of lesser length will be elided.')

    parser.add_argument(
        '--maxSequenceLen', type=int, default=None,
        help='sequences of greater length will be elided.')

    parser.add_argument(
        '--taxonomy', default=None,
        help='a string of the taxonomic group on which should be '
        'filtered. eg "Vira" will filter on viruses.')

    # Args for filtering on TitlesAlignments.
    parser.add_argument(
        '--minMatchingReads', type=int, default=None,
        help='sequences that are matched by fewer reads will be elided.')

    parser.add_argument(
        '--maxMatchingReads', type=int, default=None,
        help='sequences that are matched by more reads will be elided.')

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

    args = parser.parse_args()
    reads = FastaReads(args.fasta)
    readsAlignments = BlastReadsAlignments(reads, args.json)

    # Convert white/blacklists lists to sets.
    if args.whitelist is not None:
        args.whitelist = set(args.whitelist)
    if args.blacklist is not None:
        args.blacklist = set(args.blacklist)

    readsAlignments.filter(
        minSequenceLen=args.minSequenceLen,
        maxSequenceLen=args.maxSequenceLen,
        minStart=args.minStart,
        maxStop=args.maxStop,
        oneAlignmentPerRead=args.oneAlignmentPerRead,
        maxHspsPerHit=args.maxHspsPerHit,
        scoreCutoff=args.scoreCutoff,
        whitelist=args.whitelist,
        blacklist=args.blacklist,
        titleRegex=args.titleRegex,
        negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter,
        taxonomy=args.taxonomy)

    reads = Reads()
    count = 0

    if (args.minMatchingReads is None and
            args.maxMatchingReads is None and args.minMedianScore is None and
            args.withScoreBetterThan is None and args.minNewReads is None):
        # No need to collect into titles, just get the read ids from
        # the matching alignments.
        for readAlignment in readsAlignments:
            reads.add(readAlignment.read)
            count += 1
    else:
        # We need to collect alignments into titles.
        titlesAlignments = TitlesAlignments(readsAlignments).filter(
            minMatchingReads=args.minMatchingReads,
            maxMatchingReads=args.maxMatchingReads,
            minMedianScore=args.minMedianScore,
            withScoreBetterThan=args.withScoreBetterThan,
            minNewReads=args.minNewReads)

        for titleAlignments in titlesAlignments.values():
            for alignment in titleAlignments.alignments:
                reads.add(alignment.read)
                count += 1

    reads.save(sys.stdout)
    print('Found %d matching reads.' % count, file=sys.stderr)

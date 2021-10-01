#!/usr/bin/env python

import sys
import argparse
from json import dump

from dark.sam import DistanceMatrix


def main(args):
    """
    Print a JSON object with reference to read scores.

    @param args: A C{Namespace} object returned by argparse, containing
        argument values passed on the command line.
    """
    dm = DistanceMatrix()
    for filename in args.samFile:
        dm.addFile(filename, scoreTag=args.scoreTag)

    if args.minMatchingReads is None:
        wantedReferenceIds = sorted(set(dm.scores))
    else:
        wantedReferenceIds = sorted(set(
            referenceId for (referenceId, reads) in dm.scores.items()
            if len(reads) >= args.minMatchingReads))
        if args.verbose:
            print(f'Found {len(wantedReferenceIds)} references with at least '
                  f'{args.minMatchingReads} matching reads.', file=sys.stderr)

    if args.verbose:
        for referenceId in wantedReferenceIds:
            nReads = len(dm.scores[referenceId])
            print(f'Reference {referenceId!r} matched {nReads} reads.')

    dump(dict((id_, dm.scores[id_]) for id_ in wantedReferenceIds),
         sys.stdout, sort_keys=True, indent=4)
    print()


def makeParser():
    """
    Make a command-line argument parser.

    @return: An C{argparse.ArgumentParser} instance.
    """
    parser = argparse.ArgumentParser(
        description=('Print a JSON object containing reference to read '
                     'distances extracted from a SAM file.'))

    parser.add_argument(
        '--samFile', action='append', required=True,
        help='The SAM file(s) to load. May be repeated.')

    parser.add_argument(
        '--minMatchingReads', type=int,
        help=('The minimum number of reads that must match a reference for it '
              'to be included.'))

    parser.add_argument(
        '--scoreTag',
        help=('The score tag to use for the alignment score. If not given, '
              '1 will be used to indicate that a read matched a reference '
              '(non-matches are not included). The default is no score tag, '
              'which is not that useful. A good choice is "AS", for the '
              'alignment score, but that has to be present in the SAM file, '
              'which means that the aligner (bowtie2, bwa, etc. has to have '
              'produced such a tag.'))

    parser.add_argument(
        '--verbose', action='store_true',
        help='Print extra information.')

    return parser


if __name__ == '__main__':
    main(makeParser().parse_args())

#!/usr/bin/env python

from __future__ import print_function

import sys

from dark.fasta import FastaReads


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin and a set of filtering criteria ',
                     'write filtered FASTA to stdout.'))

    parser.add_argument(
        '--minLength', type=int,
        help='The minimum sequence length')

    parser.add_argument(
        '--maxLength', type=int,
        help='The maximum sequence length')

    parser.add_argument(
        '--removeGaps', action='store_true', default=False,
        help=("If True, gap ('-') characters in sequences will be removed."))

    parser.add_argument(
        '--whitelist', action='append', default=None,
        help='sequence titles (ids) that should be whitelisted')

    parser.add_argument(
        '--blacklist', action='append', default=None,
        help='sequence titles (ids) that should be blacklisted')

    parser.add_argument(
        '--titleRegex', default=None,
        help='a regex that sequence titles (ids) must match.')

    parser.add_argument(
        '--negativeTitleRegex', default=None,
        help='a regex that sequence titles (ids) must not match.')

    parser.add_argument(
        '--truncateTitlesAfter', default=None,
        help=('a string that sequence titles (ids) will be truncated beyond. '
              'If the truncated version of a title has already been seen, '
              'that title will be skipped.'))

    parser.add_argument(
        '--indices', type=int, action='append', default=None,
        help='sequence indices that should be returned (zero based)')

    parser.add_argument(
        '--head', type=int, default=None, metavar='N',
        help='only the first N sequences will be printed.')

    parser.add_argument(
        '--removeDuplicates', action='store_true', default=False,
        help=('If True, duplicate sequences will be removed. The first '
              'occurrence is kept.'))

    args = parser.parse_args()
    kept = 0
    reads = FastaReads(sys.stdin, checkAlphabet=False)

    for seq in reads.filter(
            minLength=args.minLength,
            maxLength=args.maxLength,
            removeGaps=args.removeGaps,
            whitelist=set(args.whitelist) if args.whitelist else None,
            blacklist=set(args.blacklist) if args.blacklist else None,
            titleRegex=args.titleRegex,
            negativeTitleRegex=args.negativeTitleRegex,
            truncateTitlesAfter=args.truncateTitlesAfter,
            indices=set(args.indices) if args.indices else None,
            head=args.head, removeDuplicates=args.removeDuplicates):
        kept += 1
        print(seq.toString('fasta'), end=' ')

    print('Read %d sequences, kept %d.' % (len(reads), kept), file=sys.stderr)

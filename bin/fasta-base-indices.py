#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import argparse

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions, unambiguousBases)

parser = argparse.ArgumentParser(
    description=(
        'Given FASTA on standard input, and bases to look for, write the '
        '1-based indices of where the bases occur in all sequences to '
        'standard output. If standard output is a terminal, these will be '
        'space separated, else newline separated. Use --any to print '
        'indices that match any sequence.'),
    epilog=(
        'This can be used to find all columns of a FASTA multiple '
        'sequence alignment that contain gaps, ambiguous nucleotides or '
        'AAs, or to find all columns that do not contain such things. '
        'Note that if sequences are of uneven lengths and --any is not '
        'specified, only indices up to the length of the shortest input '
        'sequence can be printed (i.e., there is a strict interpretation '
        'of all sequences needing to have a matching base at an index: '
        'if the index does not exist in even one sequence, then no base '
        'can occur in all sequences at that index).'))

parser.add_argument(
    '--bases',
    help=('The sequence bases whose indices should be printed. If not '
          'specified, this will be the defined set of bases for the input '
          'sequence type (i.e., "ACGT" for DNA). This will have the effect of '
          'printing the indices for which any sequence has an ambiguous or '
          'missing base.'))

parser.add_argument(
    '--matchCase', default=False, action='store_true',
    help='If specified, sequence case will be considered in matching.')

parser.add_argument(
    '--any', default=False, action='store_true',
    help=('If specified, print indices of bases that match any sequence, '
          'otherwise (the default) indices are only printed if they match '
          'all sequences.'))

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)

if args.bases is None:
    # No target bases were given. Use the set of unambiguous bases for
    # the read type. Unless --any has been given, this will print
    # indices in which no ambiguous bases or gaps appear in any
    # sequence.
    targets = unambiguousBases[args.readClass]
else:
    targets = set(args.bases)

indices = reads.indicesMatching(targets, args.matchCase, args.any)
nIndices = len(indices)

if nIndices:
    separator = ' ' if os.isatty(1) else '\n'
    # Add one to indices to maximize happy humans.
    print(separator.join(map(lambda x: str(x + 1), sorted(indices))))

print('Found %d %s where %s a base from the set {%s}.' %
      (nIndices,
       'index' if nIndices == 1 else 'indices',
       'any sequence has' if args.any else 'all sequences have',
       ', '.join(sorted(targets))),
      file=sys.stderr)

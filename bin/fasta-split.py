#!/usr/bin/env python

import sys
import argparse
from os.path import exists
from os import mkdir, rename
from math import log10
from pathlib import Path

from dark.reads import (addFASTACommandLineOptions,
                        parseFASTACommandLineOptions, Reads)
from dark.utils import take


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        'Split sequences in a FASTA/Q file into numerically named files, '
        'each containing a given number of sequences.'))

parser.add_argument(
    '--outDir', default='.',
    help='The directory to make the files in.')

parser.add_argument(
    '--verbose', action='store_true',
    help='If given, print sequence ids as they are processed.')

parser.add_argument(
    '--count', type=int, default=100,
    help='The number of sequences to put in each file.')

parser.add_argument(
    '--noLeadingZeroes', action='store_true',
    help='If given, numeric filenames will not have leading zeroes.')

parser.add_argument(
    '--force', action='store_true',
    help='If given, overwrite pre-existing files.')

parser.add_argument(
    '--saveAs', choices=('fasta', 'fastq', 'fasta-ss'),
    help=('The output format. The default is to match the input format, '
          'so there is usually no need to specify this option. It can be '
          'used to force conversion from FASTQ to FASTA'))

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)

if not exists(args.outDir):
    mkdir(args.outDir)

saveAs = (
    args.saveAs or
    (args.fasta and 'fasta') or
    (args.fastq and 'fastq') or
    (args.fasta_ss and 'fasta-ss'))

# Note: we may be reading the FASTA input from stdin, so we cannot read it
# more than once (and I don't want to store it all because it may be very
# large). That's why we do a second phase of processing to renumber the
# files we created (if --noLeadingZeroes is not used).

outDir = Path(args.outDir)

count = 0
for count, reads in enumerate(
        take(parseFASTACommandLineOptions(args), args.count), start=1):
    filename = outDir / f'{count}.fasta'
    if not args.force and filename.exists():
        print(f'Will not overwrite pre-existing file {str(filename)!r}. '
              f'Use --force to make me. Exiting.', file=sys.stderr)
        sys.exit(1)
    if args.verbose:
        print(f'Writing {filename}')
    with open(filename, 'w') as fp:
        Reads(reads).save(fp, saveAs)

# Rename numeric filenames to have leading zeroes.
if count and not args.noLeadingZeroes:
    width = int(log10(count)) + 1
    if width > 1:
        if args.verbose:
            print('Renaming to add leading zeroes.')
        for i in range(1, count + 1):
            old = str(i)
            new = '%0*d' % (width, i)
            if old == new:
                # We've reached the point where the new names have as many
                # digits as the old, so we can stop.
                break
            else:
                oldFilename = outDir / f'{old}.fasta'
                newFilename = outDir / f'{new}.fasta'
                if newFilename.exists() and not args.force:
                    print(f'Will not overwrite pre-existing file '
                          f'{str(newFilename)!r}. Use --force to make me. '
                          f'Exiting.', file=sys.stderr)
                    sys.exit(1)
                if args.verbose:
                    print(f'  {str(oldFilename)} -> {str(newFilename)}')
                rename(oldFilename, newFilename)

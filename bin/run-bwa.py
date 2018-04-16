#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from os.path import exists, basename

from dark.process import Executor

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Run bwa on a FASTA file. Optionally convert the result to '
                 'BAM, sorting, and indexing.'))

parser.add_argument(
    '--bwaIndex', required=True,
    help='The name of the BWA index (created with bwa index).')

parser.add_argument(
    '--bwaArgs', default='mem',
    help='All arguments to be passed to BWA.')

parser.add_argument(
    '--fastaFile', required=True,
    help=('The name of the FASTA file whose sequences should be aligned with '
          'the BWA index.'))

parser.add_argument(
    '--base',
    help=('The base name of files to create. Suffixes such as .sam and .bam '
          'will be added. If not given, the basename of the --fastaFile name '
          'will be used, stripped of its final suffix.'))

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help='If specified, show the commands that were (or would be) executed.')

parser.add_argument(
    '--noBAM', default=False, action='store_true',
    help='If specified, do not convert SAM to BAM.')

parser.add_argument(
    '--noSort', default=False, action='store_true',
    help='If specified, do not sort the BAM.')

parser.add_argument(
    '--noIndex', default=False, action='store_true',
    help='If specified, do not index the sorted BAM.')

parser.add_argument(
    '--noClean', default=False, action='store_true',
    help='If specified, do not remove intermediate .sam and .bam files.')

parser.add_argument(
    '--force', default=False, action='store_true',
    help='If specified, overwrite pre-existing output files.')

parser.add_argument(
    '--dryRun', default=False, action='store_true',
    help='If specified, do not run commands, just print what would be done.')

args = parser.parse_args()

if args.base is None:
    fields = basename(args.fastaFile).rsplit('.', 1)
    if len(fields) < 2:
        print('No --base argument was given and the --fastaFile argument '
              'does not have a .suffix that can be stripped.',
              file=sys.stderr)
        sys.exit(1)
    base = fields[0]
else:
    base = args.base


samFile = base + '.sam'
bamFile = base + '.bam'
sortedBamFile = base + '-sorted.bam'

if not (args.force or args.dryRun):
    existing = []
    for filename in samFile, bamFile, sortedBamFile:
        if exists(filename):
            existing.append(filename)
    if existing:
        print('Will not overwrite pre-existing file%s %s. '
              'Use --force to make me.' % (
                  '' if len(existing) == 1 else 's',
                  ', '.join(existing)),
              file=sys.stderr)
        sys.exit(2)

e = Executor(args.dryRun)

e.execute("bwa %s '%s' '%s' > '%s'" % (
    args.bwaArgs, args.bwaIndex, args.fastaFile, samFile))

if not args.noBAM:
    e.execute("samtools view -b < '%s' > '%s'" % (samFile, bamFile))

    if not args.noClean:
        e.execute("rm '%s'" % samFile)

    if not args.noSort:
        e.execute("samtools sort '%s' > '%s'" % (bamFile, sortedBamFile))

        if not args.noClean:
            e.execute("rm '%s'" % bamFile)

        if not args.noIndex:
            e.execute("samtools index '%s'" % sortedBamFile)

if args.dryRun or args.verbose:
    print('\n'.join(e.log))

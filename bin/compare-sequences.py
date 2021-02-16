#!/usr/bin/env python

from __future__ import print_function, division

import sys
import argparse
from math import log10
import multiprocessing

from dark.aligners import mafft, needle
from dark.dna import compareDNAReads, matchToString, AMBIGUOUS
from dark.reads import (Reads, addFASTACommandLineOptions,
                        parseFASTACommandLineOptions)
from dark.utils import parseRangeExpression

MAFFT_DEFAULT_ARGS = '--globalpair --maxiterate 1000 --preservecase'
MAFFT_ALGORITHMS_URL = (
    'https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html')
NEEDLE_DEFAULT_ARGS = 'auto'


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Compare two sequences.')

parser.add_argument(
    '--index1', type=int, default=1,
    help='The (1-based) index in the input of the first sequence.')

parser.add_argument(
    '--index2', type=int, default=2,
    help='The (1-based) index in the input of the second sequence.')

parser.add_argument(
    '--align', default=False, action='store_true',
    help=('If given, use mafft (the default) or needle (according to the '
          'algorithm selected by --aligner) to align the two sequences.'))

parser.add_argument(
    '--aligner', default='mafft', choices=('mafft', 'needle'),
    help='The alignment algorithm to use.')

parser.add_argument(
    '--alignerOptions',
    help=('Optional arguments to pass to the alignment algorithm. If the '
          'aligner is mafft, the default options are %r. If needle, "%s". '
          'Do not try to set the number of threads here - use the --threads '
          'argument instead. If you are using mafft, see %s for some possible '
          'option combinations.' %
          (MAFFT_DEFAULT_ARGS, NEEDLE_DEFAULT_ARGS, MAFFT_ALGORITHMS_URL)))

parser.add_argument(
    '--threads', type=int, default=multiprocessing.cpu_count(),
    help=('The number of threads to use when running the aligner (if --align '
          'is used and the alignment algorithm can make use of multiple '
          'threads (mafft can, needle cannot)).'))

parser.add_argument(
    '--alignmentFile',
    help='The file to save the alignment to (implies --align).')

parser.add_argument(
    '--strict', default=False, action='store_true',
    help='If given, do not allow ambiguous nucleotide symbols to match.')

parser.add_argument(
    '--quiet', dest='verbose', default=True, action='store_false',
    help=('Do not print information about aligning, or falling back to '
          'stretcher.'))

parser.add_argument(
    '--noGapLocations', dest='includeGapLocations', action='store_false',
    default=True,
    help='Do not indicate the (1-based) locations of sequence gaps.')

parser.add_argument(
    '--sites',
    help=('Specify (1-based) sequence sites to keep. All other sites will '
          'be ignored. The sites must be given in the form e.g., '
          '24,100-200,260.'))

parser.add_argument(
    '--showDiffs', default=False, action='store_true',
    help='Print (1-based) sites where the sequence nucleotides differ.')

parser.add_argument(
    '--showAmbiguous', default=False, action='store_true',
    help=('Print (1-based) sites where either sequence has an ambiguous '
          'nucleotide code.'))

addFASTACommandLineOptions(parser)
args = parser.parse_args()

keepSequences = set([args.index1 - 1, args.index2 - 1])

reads = list(parseFASTACommandLineOptions(args).filter(
    keepSequences=keepSequences))

if len(reads) == 1:
    if len(keepSequences) == 1:
        # This is ok, they want to compare a sequence with itself.
        reads = Reads([reads[0], reads[0]])
    else:
        print('Could not find both requested sequence indices. Exiting.')
        sys.exit(1)
elif len(reads) != 2:
    print('Could not find both requested sequence indices. Exiting.')
    sys.exit(1)

if args.alignmentFile:
    args.align = True

if args.align:
    len1, len2 = map(len, reads)
    if len1 == len2:
        print('Pre-alignment, sequence lengths were identical: %d' % len1)
    else:
        print('Pre-alignment, sequence lengths: %d, %d (difference %d)' % (
            len1, len2, abs(len1 - len2)))

    if args.aligner == 'mafft':
        # Be careful in examining args.alignerOptions because we want the
        # user to be able to pass an empty string (so check against None
        # before deciding to use the default.)
        options = (MAFFT_DEFAULT_ARGS if args.alignerOptions is None
                   else args.alignerOptions)
        reads = mafft(reads, args.verbose, options=options,
                      threads=args.threads)
    else:
        assert args.aligner == 'needle'
        # Be careful in examining args.alignerOptions because we want the
        # user to be able to pass an empty string (so check against None
        # before deciding to use the default.)
        options = (NEEDLE_DEFAULT_ARGS if args.alignerOptions is None
                   else args.alignerOptions)
        reads = needle(reads, args.verbose, options=options)

    if args.alignmentFile:
        assert reads.save(args.alignmentFile) == 2

offsets = (parseRangeExpression(args.sites, convertToZeroBased=True)
           if args.sites else None)

read1, read2 = reads
len1, len2 = map(len, reads)
identicalLengths = len1 == len2

# Sanity check.
if args.align:
    assert identicalLengths

match = compareDNAReads(read1, read2, matchAmbiguous=(not args.strict),
                        offsets=offsets)

x = 'Post-alignment, sequence' if args.align else 'Sequence'
if identicalLengths:
    print('%s lengths are identical: %s' % (x, len1))
else:
    print('%s lengths: %d, %d (difference %d)' % (x, len1, len2,
                                                  abs(len1 - len2)))

print(matchToString(match, read1, read2, matchAmbiguous=(not args.strict),
                    offsets=offsets,
                    includeGapLocations=args.includeGapLocations))

if args.showDiffs:
    # Print all sites where the sequences differ.
    width = int(log10(max(len1, len2))) + 1
    headerPrinted = False
    for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence),
                                  start=1):
        if a != b:
            if not headerPrinted:
                print('Differences (site, %s, %s):' % (read1.id, read2.id))
                headerPrinted = True
            print('  %*d %s %s' % (width, site, a, b))

    if not headerPrinted:
        print('No sequence differences found.')

if args.showAmbiguous:
    width = int(log10(max(len1, len2))) + 1
    headerPrinted = False
    for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence),
                                  start=1):
        if len(AMBIGUOUS.get(a, '')) > 1 or len(AMBIGUOUS.get(b, '')) > 1:
            if not headerPrinted:
                print('Ambiguities (site, %s, %s):' % (read1.id, read2.id))
                headerPrinted = True
            print('  %*d %s %s' % (width, site, a, b))

    if not headerPrinted:
        print('No sequence ambiguities found.')

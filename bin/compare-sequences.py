#!/usr/bin/env python

from __future__ import print_function, division

import sys
import argparse
from os.path import join
from math import log10

from dark.dna import compareDNAReads, matchToString
from dark.fasta import FastaReads
from dark.reads import (Reads, addFASTACommandLineOptions,
                        parseFASTACommandLineOptions)
from dark.process import Executor
from dark.utils import parseRangeString


def needle(reads):
    """
    Run a Needleman-Wunsch alignment and return the two sequences.

    @param reads: An iterable of two reads.
    @return: A C{Reads} instance with the two aligned sequences.
    """
    from tempfile import mkdtemp
    from shutil import rmtree

    dir = mkdtemp()

    file1 = join(dir, 'file1.fasta')
    with open(file1, 'w') as fp:
        print(reads[0].toString('fasta'), end='', file=fp)

    file2 = join(dir, 'file2.fasta')
    with open(file2, 'w') as fp:
        print(reads[1].toString('fasta'), end='', file=fp)

    out = join(dir, 'result.fasta')

    Executor().execute("needle -asequence '%s' -bsequence '%s' -auto "
                       "-outfile '%s' -aformat fasta" % (
                           file1, file2, out))

    # Use 'list' in the following to force reading the FASTA from disk.
    result = Reads(list(FastaReads(out)))
    rmtree(dir)

    return result


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Compare two sequences.'))

parser.add_argument(
    '--index1', type=int, default=1,
    help='The (1-based) index in the input of the first sequence.')

parser.add_argument(
    '--index2', type=int, default=2,
    help='The (1-based) index in the input of the second sequence.')

parser.add_argument(
    '--align', default=False, action='store_true',
    help='If given, use needle to do a global alignment of the two sequences.')

parser.add_argument(
    '--alignmentFile',
    help='The file to save the alignment to (implies --align).')

parser.add_argument(
    '--strict', default=False, action='store_true',
    help='If given, do not allow ambiguous nucleotide symbols to match')

parser.add_argument(
    '--sites',
    help=('Specify (1-based) sequence sites to keep. All other sites will '
          'be ignored. The sites must be given in the form e.g., '
          '24,100-200,260.'))

parser.add_argument(
    '--showDiffs', default=False, action='store_true',
    help='Print (1-based) sites where the sequence nucleotides differ.')

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
        print('Pre-alignment, sequence lengths were identical: %s' % len1)
    else:
        print('Pre-alignment, sequence lengths: %d, %d (difference %d)' % (
            len1, len2, abs(len1 - len2)))

    # Align.
    reads = needle(reads)

    if args.alignmentFile:
        assert reads.save(args.alignmentFile) == 2

offsets = (parseRangeString(args.sites, convertToZeroBased=True)
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
                    offsets=offsets))

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

#!/usr/bin/env python

from __future__ import print_function, division

import sys
import argparse
from os.path import join

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.reads import (Reads, addFASTACommandLineOptions,
                        parseFASTACommandLineOptions)
from dark.subprocess import Executor


def needle(reads):
    """
    Run a Needleman-Wunsch alignment and return the two sequences.

    @param reads: An interable of two reads.
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


def pct(a, b):
    """
    What percent of a is b?

    @param a: a numeric value.
    @param b: a numeric value.
    @return: the C{float} percentage.
    """
    return 100.0 * a / b if b else 0.0


def pp(mesg, count, len1, len2=None):
    """
    Print a message followed by an integer count and a percentage (or
    two, if the sequence lengths are unequal).

    @param mesg: a C{str} message.
    @param count: a numeric value.
    @param len1: the C{int} length of sequence 1.
    @param len2: the C{int} length of sequence 2. If not given, will
        default to C{len1}.
    """
    if count == 0:
        print('%s: %d' % (mesg, count))
    else:
        len2 = len2 or len1
        if len1 == len2:
            print('%s: %d (%.2f%%)' % (mesg, count, pct(count, len1)))
        else:
            print('%s: %d (%.2f%% of sequence 1; %.2f%% of sequence 2)' % (
                mesg, count, pct(count, len1), pct(count, len2)))


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

addFASTACommandLineOptions(parser)
args = parser.parse_args()

indices = set([args.index1 - 1, args.index2 - 1])

reads = list(parseFASTACommandLineOptions(args).filter(indices=indices))

if len(reads) == 1:
    if len(indices) == 1:
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

read1, read2 = reads
len1, len2 = map(len, reads)
identicalLengths = len1 == len2

# Sanity check.
if args.align:
    assert identicalLengths

result = compareDNAReads(read1, read2, matchAmbiguous=(not args.strict))

match = result['match']
identicalMatchCount = match['identicalMatchCount']
ambiguousMatchCount = match['ambiguousMatchCount']
gapMismatchCount = match['gapMismatchCount']
gapGapMismatchCount = match['gapGapMismatchCount']
nonGapMismatchCount = match['nonGapMismatchCount']

x = 'Post-alignment, sequence' if args.align else 'Sequence'
if identicalLengths:
    print('%s lengths are identical: %s' % (x, len1))
else:
    print('%s lengths: %d, %d (difference %d)' % (x, len1, len2,
                                                  abs(len1 - len2)))

pp('Exact matches', identicalMatchCount, len1, len2)
pp('Ambiguous matches', ambiguousMatchCount, len1, len2)
if ambiguousMatchCount and identicalMatchCount:
    anyMatchCount = identicalMatchCount + ambiguousMatchCount
    pp('Exact or ambiguous matches', anyMatchCount, len1, len2)
mismatchCount = gapMismatchCount + gapGapMismatchCount + nonGapMismatchCount
pp('Mismatches', mismatchCount, len1, len2)
conflicts = 'conflicts or ambiguities' if args.strict else 'conflicts'
pp('  Not involving gaps (i.e., %s)' % conflicts,
   nonGapMismatchCount, len1, len2)
pp('  Involving a gap in one sequence', gapMismatchCount, len1, len2)
pp('  Involving a gap in both sequences', gapGapMismatchCount, len1, len2)

for read, index, key in zip(reads, (args.index1, args.index2),
                            ('read1', 'read2')):
    print('Sequence index %d:' % index)
    print('  Id: %s' % read.id)
    length = len(read)
    print('  Length: %d' % length)
    gapCount = len(result[key]['gapOffsets'])
    pp('  Gaps', gapCount, length)
    ambiguousCount = len(result[key]['ambiguousOffsets'])
    pp('  Ambiguous', ambiguousCount, length)
    extraCount = result[key]['extraCount']
    if extraCount:
        pp('  Extra nucleotides at end', extraCount, length)

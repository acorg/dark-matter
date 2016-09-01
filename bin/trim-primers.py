#!/usr/bin/env python

from __future__ import print_function

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from dark.sequence import findPrimerBidiLimits


def trimPrimers(primer, verbose):
    """
    @param primer: A BioPython C{Bio.Seq} primer sequence.
    @param verbose: A C{bool}, if C{True} output additional information about
        how often and where primers were found.
    """
    reads = []
    absentCount = forwardCount = reverseCount = count = 0
    for seqRecord in SeqIO.parse(sys.stdin, 'fasta'):
        count += 1
        start, end = findPrimerBidiLimits(primer, seqRecord.seq)
        if start == 0:
            if end == len(seqRecord):
                absentCount += 1
            else:
                reverseCount += 1
        else:
            forwardCount += 1
            if end != len(seqRecord):
                reverseCount += 1
        reads.append(seqRecord[start:end])

    if verbose:
        print((
            'Read %d sequences. Found forward: %d, '
            'Found reversed: %d, Absent: %d') % (
            count, forwardCount, reverseCount, absentCount), file=sys.stderr)

    SeqIO.write(reads, sys.stdout, 'fasta')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin, look for a primer sequence '
                     'and write trimmed FASTA (after the primer) to stdout.'))

    parser.add_argument('primer', help='the primer sequence')
    parser.add_argument('--verbose', type=bool, default=False,
                        help='If True, print information on found primers')

    args = parser.parse_args()

    trimPrimers(Seq(args.primer.upper(), IUPAC.unambiguous_dna),
                args.verbose)

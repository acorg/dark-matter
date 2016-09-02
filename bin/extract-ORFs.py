#!/usr/bin/env python

"""
Read DNA or AA FASTA from stdin and print AA ORF-only FASTA to stdout. I.e.,
the output FASTA will contain one sequence for each ORF found in each
sequence on stdin.

If a minimal ORF length is given, only print ORFs of at least that length
unless we are asked to allow ORFs that are open at at least one end.

Note that no start or stop codons will appear in the output. If you are
wanting to simply translate DNA FASTA to AA FASTA on a sequence-by-sequence
basis leaving in start/stop codons, use dna-to-aa.py instead.
"""

from __future__ import print_function

import sys
from os.path import basename
import argparse

from Bio.Data.CodonTable import TranslationError

from dark.fasta import FastaReads
from dark.reads import AARead, RNARead, DNARead

TYPE = {
    'aa': AARead,
    'dna': DNARead,
    'rna': RNARead,
}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert DNA or AA FASTA to AA ORF-only FASTA',
        epilog='Given DNA or AA FASTA on stdin, output AA ORF-only FASTA to '
        'stdout. Optionally, filter by minimum required ORF length and allow '
        'the output of short ORFs that are open.'
    )

    parser.add_argument(
        '--allowOpenORFs', default=False, action='store_true',
        help='If True, ORFs that do not meet the length requirement '
        '(as specified by --minORFLength) will be output as long as '
        'they are open.')

    parser.add_argument(
        '--minORFLength', metavar='LEN', type=int, default=None,
        help='Only ORFs of at least this length will be written to stdout.')

    parser.add_argument(
        '--type', default='dna', choices=TYPE.keys(),
        help='The type of the bases in the stdin FASTA.')

    args = parser.parse_args()
    write = sys.stdout.write
    minORFLength = args.minORFLength
    allowOpenORFs = args.allowOpenORFs

    reads = FastaReads(sys.stdin, readClass=TYPE[args.type])

    for read in reads:
        if args.type == 'aa':
            # No need to translate to AA, since our input is already AA.
            translations = [read]
        else:
            translations = read.translations()
        try:
            for translation in translations:
                for orf in translation.ORFs():
                    if minORFLength is None or len(orf) >= minORFLength or (
                            allowOpenORFs and (orf.openLeft or orf.openRight)):
                        write(orf.toString('fasta'))
        except TranslationError as error:
            print('Could not translate read %r sequence %r (%s).\nDid you '
                  'forget to run %s with "--type aa"?' %
                  (read.id, read.sequence, error, basename(sys.argv[0])),
                  file=sys.stderr)
            sys.exit(1)

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

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert DNA or AA FASTA to AA ORF-only FASTA',
        epilog=('Given DNA or AA FASTA on stdin, output AA ORF-only FASTA to '
                'stdout. Optionally, filter by minimum required ORF length '
                'and allow the output of short ORFs that are open.'))

    parser.add_argument(
        '--allowOpenORFs', default=False, action='store_true',
        help=('If True, ORFs that do not meet the length requirement '
              '(as specified by --minORFLength) will be output as long as '
              'they are open.'))

    parser.add_argument(
        '--minORFLength', metavar='LEN', type=int, default=None,
        help='Only ORFs of at least this length will be written to stdout.')

    addFASTACommandLineOptions(parser)

    args = parser.parse_args()
    allowOpenORFs = args.allowOpenORFs
    write = sys.stdout.write
    minORFLength = args.minORFLength
    reads = parseFASTACommandLineOptions(args)
    aa = args.readClass in ('AARead', 'AAReadORF', 'AAReadWithX', 'SSAARead',
                            'SSAAReadWithX', 'TranslatedRead')

    if aa:
        def translations(read):
            return (read,)
    else:
        def translations(read):
            return read.translations()

    for read in reads:
        try:
            for translation in translations(read):
                for orf in translation.ORFs():
                    if minORFLength is None or len(orf) >= minORFLength or (
                            allowOpenORFs and (orf.openLeft or orf.openRight)):
                        write(orf.toString('fasta'))
        except TranslationError as error:
            print('Could not translate read %r sequence %r (%s).' %
                  (read.id, read.sequence, error), file=sys.stderr)
            if not aa:
                print('Did you forget to run %s with "--readClass AARead"?' %
                      (basename(sys.argv[0])), file=sys.stderr)
            sys.exit(1)

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

    parser.add_argument(
        '--maxORFLength', metavar='LEN', type=int, default=None,
        help='Only ORFs of a maximum of this length will be written to'
        'stdout.')

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '--kozakInfo', default=False, action='store_true',
        help='A file with Kozak consensus information will be written'
        'containing all Kozak consensus sequences. Only applicable if DNA'
        'reads given.')

    parser.add_argument(
        '--kozakOnly', default=False, action='store_true',
        help='Only ORFs that also have a Kozak consensus will be written to'
        'stdout. A file with Kozak consensus information will be written'
        'containing only the Kozak consensus sequences that correspond to an'
        'ORF. Only applicable if DNA reads given.')

    parser.add_argument(
        '--kozakInfoFile', default='kozakInfo.out',
        help='The name of the file to which the Kozak consensus information'
        'will be written.')

    addFASTACommandLineOptions(parser)

    args = parser.parse_args()
    allowOpenORFs = args.allowOpenORFs
    minORFLength = args.minORFLength
    maxORFLength = args.maxORFLength
    kozakInfo = args.kozakInfo
    kozakInfoFile = args.kozakInfoFile
    kozakOnly = args.kozakOnly
    reads = parseFASTACommandLineOptions(args)
    aa = args.readClass in ('AARead', 'AAReadORF', 'AAReadWithX', 'SSAARead',
                            'SSAAReadWithX', 'TranslatedRead')

    if aa:
        def translations(read):
            return (read,)
    else:
        def translations(read):
            return read.translations()

    if kozakInfo or kozakOnly:
        if aa:
            # Should the code break in this case?
            print('Did you pass NA sequences? Kozak consensuses can only'
                  'be found in a nucleotide sequence.', file=sys.stderr)
        else:
            from dark.dna import findKozakConsensus

        def writeToKozakOut(kozakread, kozakInfoFile):
            """
            Writes out information about a Kozak sequence stored in kozakread
            to the file name given in kozakInfoFile.
            """
            kozakout = open(kozakInfoFile, 'a')
            kozakout.write('Read ID: ' + kozakread.id)
            kozakout.write('\n')
            kozakout.write('Kozak sequence: ' + kozakread.sequence)
            kozakout.write('\n')
            kozakout.write('Offset of the Kozak sequence A: ' +
                           str(kozakread.start + 6))
            kozakout.write('\n')
            kozakout.write('Quality of the Kozak consensus: ' +
                           str(kozakread.kozakQuality) + ' %')
            kozakout.write('\n')
            kozakout.write('\n')
            kozakout.close()

    for read in reads:
        try:
            for translation in translations(read):
                for orf in translation.ORFs(allowOpenORFs):
                    # Check the length requirements, if any.
                    if ((minORFLength is None or len(orf) >= minORFLength) and
                       (maxORFLength is None or len(orf) <= maxORFLength)):
                        # If Kozak consensus information is wanted:
                        if kozakInfo and not aa:
                            kozakreads = list(findKozakConsensus(read))
                            for kozakread in kozakreads:
                                writeToKozakOut(kozakread, kozakInfoFile)
                            sys.stdout.write(orf.toString('fasta'))
                        # If only ORFs with a Kozak consensus are wanted:
                        elif kozakOnly and not aa:
                            kozakreads = list(findKozakConsensus(read))
                            for kozakread in kozakreads:
                                if kozakread.stop in (orf.start * 3 + 1,
                                                      orf.start * 3 + 2,
                                                      orf.start * 3 + 3):
                                    writeToKozakOut(kozakread, kozakInfoFile)
                                    sys.stdout.write(orf.toString('fasta'))
                        # If no Kozak Info is wanted write ORF to stdout.
                        else:
                            sys.stdout.write(orf.toString('fasta'))

        except TranslationError as error:
            print('Could not translate read %r sequence %r (%s).' %
                  (read.id, read.sequence, error), file=sys.stderr)
            if not aa:
                print('Did you forget to run %s with "--readClass AARead"?' %
                      (basename(sys.argv[0])), file=sys.stderr)
            sys.exit(1)

#!/usr/bin/env python

"""
Extract aligned (i.e., padded) queries in FASTA format from a SAM/BAM file.
"""

from __future__ import division, print_function

import sys
import argparse

from dark.sam import (
    PaddedSAM, UnequalReferenceLengthError, UnknownReference)

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Produce aligned FASTA queries from a SAM/BAM file.')

parser.add_argument(
    'samFile', help='The SAM/BAM file to read.')

parser.add_argument(
     '--minLength', type=int, default=0,
     help='Do not output anything for queries shorter than this.')

parser.add_argument(
     '--referenceName',
     help=('The name of the reference sequence to print alignments for. '
           'This is only needed if the SAM/BAM alignment was made against '
           'multiple references *and* they are not all of the same length. '
           'If there is only one reference sequence or if all reference '
           'sequences are of the same length, there is no need to provide a '
           'reference name.'))

parser.add_argument(
     '--rcSuffix', default='',
     help=('A string to add to the end of query names that are reverse '
           'complemented. This is added before the /1, /2, etc., that are '
           'added for duplicated ids (if there are duplicates and '
           '--allowDuplicateIds is not used)'))

parser.add_argument(
     '--dropSecondary', default=False, action='store_true',
     help='If given, secondary matches will not be output.')

parser.add_argument(
     '--dropSupplementary', default=False, action='store_true',
     help='If given, supplementary matches will not be output.')

parser.add_argument(
     '--allowDuplicateIds', default=False, action='store_true',
     help=('If given, repeated query ids (due to secondary or supplemental '
           'matches) will not have /1, /2, etc. appended to their ids. So '
           'repeated ids may appear in the output FASTA.'))

parser.add_argument(
     '--keepQCFailures', default=False, action='store_true',
     help=('If given, reads that are considered quality control failures will '
           'be included in the output.'))

parser.add_argument(
     '--listReferenceNames', default=False, action='store_true',
     help=('If given, the names of all reference sequences will be '
           'printed to standard output and the program will exit.'))

parser.add_argument(
     '--rcNeeded', default=False, action='store_true',
     help=('If given, queries that are flagged as matching when reverse '
           'complemented will be reverse complemented in the output. This '
           'must be used if the program that created the SAM/BAM input '
           'flags reversed matches but does not also store the reverse '
           'complemented query. The bwa program (mem and aln followed by '
           'samse) stores the queries reversed complemented if the match '
           'was, so this option is not needed for bwa. If in doubt, test the '
           'output of your matching program as this is very important!'))


args = parser.parse_args()

paddedSAM = PaddedSAM(args.samFile)

try:
    if args.listReferenceNames:
        print(paddedSAM.referencesToStr(), file=sys.stdout)
    else:
        try:
            for read in paddedSAM.queries(
                    referenceName=args.referenceName,
                    minLength=args.minLength,
                    rcSuffix=args.rcSuffix,
                    dropSecondary=args.dropSecondary,
                    dropSupplementary=args.dropSupplementary,
                    allowDuplicateIds=args.allowDuplicateIds,
                    keepQCFailures=args.keepQCFailures,
                    rcNeeded=args.rcNeeded):
                print(read.toString('fasta'), end='')
        except UnequalReferenceLengthError as e:
            raise ValueError(
                str(e) + ' So it is not clear how long the padded output '
                'FASTA sequences should be. Use --referenceName to specify '
                'which reference sequence is the one whose aligned reads you '
                'want printed. Use --listReferenceNames to see a list of '
                'reference sequence names and lengths.')
        except UnknownReference as e:
            raise ValueError(
                str(e) + ' Use --listReferenceNames to see a list of '
                'reference sequence names.')
finally:
    paddedSAM.close()

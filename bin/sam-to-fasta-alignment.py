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
     '--dropDuplicates', default=False, action='store_true',
     help=('If given, matches flagged as optical or PCR duplicates will '
           'not be output.'))

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
     '--listReferenceInsertions', default=False, action='store_true',
     help=('If given, information about reference sequence insertions will be '
           'printed to standard error. These correspond to "I" CIGAR '
           'operations that for the match would require inserting query bases '
           'into the reference. Because we cannot change the reference (in '
           'fact we typically do not have the reference in the SAM/BAM file), '
           'we cut the inserted bases out of the aligned query and save the '
           'information about what would have been inserted and where. That '
           'information is printed by this option. The output gives the '
           '0-based offset where the inserted base would be placed, followed '
           'by a list of the nucleotides that were suggested as being '
           'inserted and the number of times each nucleotide was suggested. '
           'So for example the output might contain "27: T:3, G:10" which '
           'indicates that 13 query (3 with T and 10 with G) matches would '
           'insert a nucleotide into the reference at offset 27.'))

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


def baseCountsToStr(counts):
    """
    Convert base counts to a string.

    @param counts: A C{Counter} instance.
    @return: A C{str} representation of nucleotide counts at an offset.
    """
    return ' '.join([
        ('%s:%d' % (base, counts[base])) for base in sorted(counts)])


def nucleotidesToStr(nucleotides, prefix=''):
    """
    Convert offsets and base counts to a string.

    @param nucleotides: A C{defaultdict(Counter)} instance, keyed
        by C{int} offset, with nucleotides keying the Counters.
    @param prefix: A C{str} to put at the start of each line.
    @return: A C{str} representation of the offsets and nucleotide
        counts for each.
    """
    result = []
    for offset in sorted(nucleotides):
        result.append(
            '%s%d: %s' % (prefix, offset,
                          baseCountsToStr(nucleotides[offset])))
    return '\n'.join(result)


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
                    dropDuplicates=args.dropDuplicates,
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

        if args.listReferenceInsertions:
            if paddedSAM.referenceInsertions:
                print('(0-based) insertions into the reference:\n%s' %
                      nucleotidesToStr(paddedSAM.referenceInsertions, '  '),
                      file=sys.stderr)
            else:
                print('No matches required an insertion into the reference.',
                      file=sys.stderr)
finally:
    paddedSAM.close()

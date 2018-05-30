#!/usr/bin/env python

"""
Convert SAM/BAM to an alignment in FASTA or Phylip format.
"""

from __future__ import division, print_function

import sys
import argparse
from collections import Counter

from dark.reads import Read, DNARead
from pysam import (
    AlignmentFile, CMATCH, CINS, CDEL, CREF_SKIP, CSOFT_CLIP, CHARD_CLIP, CPAD,
    CEQUAL, CDIFF)

# See https://samtools.github.io/hts-specs/SAMv1.pdf for the True/False values
# in the following and http://pysam.readthedocs.io/en/latest/api.html for the
# pysam documentation.

CONSUMES_QUERY = {
    int(CMATCH): True,
    int(CINS): True,
    int(CDEL): False,
    int(CREF_SKIP): False,
    int(CSOFT_CLIP): True,
    int(CHARD_CLIP): False,
    int(CPAD): False,
    int(CEQUAL): True,
    int(CDIFF): True,
}

CONSUMES_REFERENCE = {
    int(CMATCH): True,
    int(CINS): False,
    int(CDEL): True,
    int(CREF_SKIP): True,
    int(CSOFT_CLIP): False,
    int(CHARD_CLIP): False,
    int(CPAD): False,
    int(CEQUAL): True,
    int(CDIFF): True,
}

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Convert a SAM/BAM file to aligned FASTA.')

parser.add_argument(
    'samFile', help='The SAM/BAM file to convert.')

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
     '--dropSecondaries', default=False, action='store_true',
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
     '--preserveComplementarity', default=False, action='store_true',
     help=('If given, reads will be printed as they appear in the SAM/BAM '
           'input, regardless of whether their reverse compliment was '
           'what matched a reference. This preserves the orientation of '
           'the read as it was put into the input file, which may be '
           'desirable, but those reads (which had to be reverse complemented '
           'to match) in the output FASTA will not align with the reference.'))


def printReferences(samfile, indent=0, fp=sys.stdout):
    """
    Print a list of known reference names and their lengths.

    @param samfile: A C{pysam.AlignmentFile} instance.
    @param indent: An C{int} specifying how many spaces to indent each line.
    @param fp: The file pointer to print to.
    """
    indent = ' ' * indent
    for i in range(samfile.nreferences):
        print('%s%s (length %d)' % (
            indent, samfile.get_reference_name(i), samfile.lengths[i]),
            file=fp)


args = parser.parse_args()
samfile = AlignmentFile(args.samFile)

if args.listReferenceNames:
    printReferences(samfile)
    sys.exit(0)

if args.referenceName:
    referenceId = samfile.get_tid(args.referenceName)
    if referenceId == -1:
        print('Reference %r is not present in the SAM/BAM file. The '
              'reference sequence names are present:' %
              args.referenceName, file=sys.stderr)
        printReferences(samfile, 2)
        sys.exit(1)
    referenceLength = samfile.lengths[referenceId]
else:
    # No reference name was given. Make sure that all references have the same
    # length.
    if len(set(samfile.lengths)) != 1:
        print('Your SAM/BAM file has %d reference sequences, and their '
              'lengths (%s) are not all identical, so it is not clear how '
              'long the output FASTA sequences should be. Use '
              '--referenceName to specify which reference sequence is the '
              'one whose aligned reads you want printed.' % (
                  samfile.nreferences, ', '.join(sorted(samfile.lengths))),
              file=sys.stderr)
        printReferences(samfile, indent=2, fp=sys.stderr)
        sys.exit(1)
    referenceId = None
    referenceLength = samfile.lengths[0]

minLength = args.minLength
rcSuffix = args.rcSuffix
keepQCFailures = args.keepQCFailures
dropSecondaries = args.dropSecondaries
dropSupplementary = args.dropSupplementary
allowDuplicateIds = args.allowDuplicateIds

# Hold the count for each id so we can add /1, /2 etc to duplicate ids (unless
# --allowDuplicateIds was given).
idCount = Counter()

for read in samfile.fetch():
    query = read.query_sequence
    queryLength = len(query)
    if (read.is_unmapped or
            queryLength < minLength or
            (read.is_secondary and dropSecondaries) or
            (read.is_supplementary and dropSupplementary) or
            (read.is_qcfail and not keepQCFailures) or
            (referenceId is not None and read.reference_id != referenceId)):
        continue

    if read.is_reverse:
        query = DNARead('id', query).reverseComplement().sequence
        if rcSuffix:
            read.query_name += rcSuffix

    index = 0
    alignedSequence = []
    for operation, length in read.cigartuples:
        operation = int(operation)
        if CONSUMES_QUERY[operation]:
            alignedSequence.append(query[index:index + length])
            index += length
        elif CONSUMES_REFERENCE[operation]:
            alignedSequence.append('-' * length)

    # Sanity check that we consumed the entire query.
    assert index == queryLength

    alignedSequence = ''.join(alignedSequence)

    if read.is_reverse and args.preserveComplementarity:
        alignedSequence = DNARead(
            'id', alignedSequence).reverseComplement().sequence

    # Put '-' gap characters before and after the aligned sequence so that it
    # is offset properly and matches the length of the reference.
    paddedSequence = (
        ('-' * read.reference_start) +
        alignedSequence +
        '-' * (referenceLength -
               (read.reference_start + len(alignedSequence))))

    if allowDuplicateIds:
        suffix = ''
    else:
        count = idCount[read.query_name]
        idCount[read.query_name] += 1
        suffix = '' if count == 0 else '/%d' % count

    print(Read(
        '%s%s' % (read.query_name, suffix),
        paddedSequence).toString('fasta'), end='')

samfile.close()

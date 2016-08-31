#!/usr/bin/env python

from __future__ import print_function

import sys

from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.fastq import FastqReads


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin and a set of filtering criteria '
                     'write filtered FASTA to stdout.'))

    parser.add_argument(
        '--readClass', default='fasta', choices=('fasta', 'fastq', 'fasta-ss'),
        help=("If specified, give the input FASTA type"))

    parser.add_argument(
        '--saveAs', default=None, choices=('fasta', 'fastq', 'fasta-ss'),
        help=('The output format. The default is to match the input format, '
              'so there is usually no need to specify this option. It can be '
              'used to force conversion from FASTQ to FASTA'))

    parser.add_argument(
        '--minLength', type=int,
        help='The minimum sequence length')

    parser.add_argument(
        '--maxLength', type=int,
        help='The maximum sequence length')

    parser.add_argument(
        '--removeGaps', action='store_true', default=False,
        help=("If True, gap ('-') characters in sequences will be removed."))

    parser.add_argument(
        '--whitelist', action='append', default=None,
        help='sequence titles (ids) that should be whitelisted')

    parser.add_argument(
        '--blacklist', action='append', default=None,
        help='sequence titles (ids) that should be blacklisted')

    parser.add_argument(
        '--titleRegex', default=None,
        help='a regex that sequence titles (ids) must match.')

    parser.add_argument(
        '--negativeTitleRegex', default=None,
        help='a regex that sequence titles (ids) must not match.')

    parser.add_argument(
        '--truncateTitlesAfter', default=None,
        help=('a string that sequence titles (ids) will be truncated beyond. '
              'If the truncated version of a title has already been seen, '
              'that title will be skipped.'))

    parser.add_argument(
        '--indices', type=int, action='append', default=None,
        help='sequence indices that should be returned (zero based)')

    parser.add_argument(
        '--head', type=int, default=None, metavar='N',
        help='only the first N sequences will be printed.')

    parser.add_argument(
        '--removeDuplicates', action='store_true', default=False,
        help=('If True, duplicate sequences will be removed. The first '
              'occurrence is kept.'))

    # See the docstring for dark.reads.Reads.filter for more detail on
    # randomSubset.
    parser.add_argument(
        '--randomSubset', type=int, default=None,
        help=('An integer giving the number of sequences that should be kept. '
              'These will be selected at random.'))

    # See the docstring for dark.reads.Reads.filter for more detail on
    # trueLength.
    parser.add_argument(
        '--trueLength', type=int, default=None,
        help=('The number of reads in the FASTA input. Only to be used with '
              'randomSubset'))

    parser.add_argument(
        '--sampleFraction', type=float, default=None,
        help=('A [0.0, 1.0] C{float} indicating a fraction of the reads that '
              'should be allowed to pass through the filter. The sample size '
              'will only be approximately the product of the sample fraction '
              'and the number of reads. The sample is taken at random.'))

    args = parser.parse_args()

    if args.readClass == 'fastq':
        # TODO: FastqReads should take a checkAlphabet argument, in the way
        # that FastaReads does.
        reads = FastqReads(sys.stdin)
    elif args.readClass == 'fasta':
        reads = FastaReads(sys.stdin, checkAlphabet=False)
    else:
        # args.readClass must be fasta-ss due to the 'choices' argument
        # passed to parser.add_argument value above.
        assert args.readClass == 'fasta-ss'
        reads = SSFastaReads(sys.stdin, checkAlphabet=False)

    saveAs = args.saveAs or args.readClass

    # Check for incompatible read/write formats. We can't write FASTQ
    # unless we have FASTQ on input (else we won't have quality
    # information), and we can't write PDB FASTA with secondary structure
    # information unless we have that on input.
    if saveAs == 'fastq' and args.readClass != 'fastq':
        raise ValueError(
            'You have specified --saveAs fastq without using --readClass '
            'fastq to indicate that the input is FASTQ. Please be explicit.')
    elif saveAs == 'fasta-ss' and args.readClass != 'fasta-ss':
        raise ValueError(
            'You have specified --saveAs fasta-ss without using --readClass '
            'fasta-ss to indicate that the input is PDB FASTA. Please be '
            'explicit.')

    kept = 0

    for seq in reads.filter(
            minLength=args.minLength,
            maxLength=args.maxLength,
            removeGaps=args.removeGaps,
            whitelist=set(args.whitelist) if args.whitelist else None,
            blacklist=set(args.blacklist) if args.blacklist else None,
            titleRegex=args.titleRegex,
            negativeTitleRegex=args.negativeTitleRegex,
            truncateTitlesAfter=args.truncateTitlesAfter,
            indices=set(args.indices) if args.indices else None,
            head=args.head, removeDuplicates=args.removeDuplicates,
            randomSubset=args.randomSubset, trueLength=args.trueLength,
            sampleFraction=args.sampleFraction):
        kept += 1
        print(seq.toString(format_=saveAs), end='')

    print('Read %d sequences, kept %d.' % (len(reads), kept), file=sys.stderr)

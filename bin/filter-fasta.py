#!/usr/bin/env python

from __future__ import print_function, division

import sys

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions
from dark.utils import parseRangeString


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Given FASTA on stdin and a set of filtering criteria '
                     'write filtered FASTA to stdout.'))

    parser.add_argument(
        '--quiet', action='store_true', default=False,
        help=('If True, do not print the final sequence summary.'))

    parser.add_argument(
        '--saveAs', choices=('fasta', 'fastq', 'fasta-ss'),
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
        '--whitelist', action='append',
        help='sequence titles (ids) that should be whitelisted')

    parser.add_argument(
        '--blacklist', action='append',
        help='sequence titles (ids) that should be blacklisted')

    parser.add_argument(
        '--whitelistFile',
        help=('the name of a file that contains sequence titles (ids) that '
              'should be whitelisted, one per line'))

    parser.add_argument(
        '--blacklistFile',
        help=('the name of a file that contains sequence titles (ids) that '
              'should be blacklisted, one per line'))

    parser.add_argument(
        '--titleRegex', help='a regex that sequence titles (ids) must match.')

    parser.add_argument(
        '--negativeTitleRegex',
        help='a regex that sequence titles (ids) must not match.')

    parser.add_argument(
        '--truncateTitlesAfter',
        help=('a string that sequence titles (ids) will be truncated beyond. '
              'If the truncated version of a title has already been seen, '
              'that title will be skipped.'))

    parser.add_argument(
        '--indices', type=int, action='append',
        help='sequence indices that should be returned (zero based)')

    parser.add_argument(
        '--head', type=int, metavar='N',
        help='only the first N sequences will be printed.')

    parser.add_argument(
        '--removeDuplicates', action='store_true', default=False,
        help=('If True, duplicate reads will be removed, based only on '
              'sequence identity. The first occurrence is kept.'))

    parser.add_argument(
        '--removeDuplicatesById', action='store_true', default=False,
        help=('If True, duplicate reads will be removed, based only on '
              'read id. The first occurrence is kept.'))

    parser.add_argument(
        '--removeDescriptions', action='store_true', default=False,
        help=('If True, read id descriptions will be removed. The '
              'description is the part of a sequence id after the '
              'first whitespace (if any).'))

    # See the docstring for dark.reads.Reads.filter for more detail on
    # randomSubset.
    parser.add_argument(
        '--randomSubset', type=int,
        help=('An integer giving the number of sequences that should be kept. '
              'These will be selected at random.'))

    # See the docstring for dark.reads.Reads.filter for more detail on
    # trueLength.
    parser.add_argument(
        '--trueLength', type=int,
        help=('The number of reads in the FASTA input. Only to be used with '
              'randomSubset'))

    parser.add_argument(
        '--sampleFraction', type=float,
        help=('A [0.0, 1.0] C{float} indicating a fraction of the reads that '
              'should be allowed to pass through the filter. The sample size '
              'will only be approximately the product of the sample fraction '
              'and the number of reads. The sample is taken at random.'))

    parser.add_argument(
        '--sequenceNumbersFile',
        help=('A file of (1-based) sequence numbers to retain. Numbers must '
              'be one per line.'))

    # A mutually exclusive group for either --keepIndices or --removeIndices.
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '--keepIndices',
        help=('Specify 1-based indices to keep. All other sequence indices '
              'will be removed. The indices must be given in the form e.g., '
              '24,100-200,260'))

    group.add_argument(
        '--keepIndicesFile',
        help=('Specify a file containing 1-based indices to keep. All other '
              'sequence indices will be removed. Lines in the file must be '
              'given in the form e.g., 24,100-200,260'))

    group.add_argument(
        '--removeIndices',
        help=('Specify 1-based indices to remove. All other sequence indices '
              'will be kept. The indices must be given in the form e.g., '
              '24,100-200,260'))

    group.add_argument(
        '--removeIndicesFile',
        help=('Specify a file containing 1-based indices to remove. All other '
              'sequence indices will be kept. Lines in the file must be given '
              'in the form e.g., 24,100-200,260'))
    parser.add_argument(
        '--checkResultCount', type=int,
        help=('The number of reads expected in the output. If this number is '
              'not seen, the script exits with status 1 and an error '
              'message is printed unless --quiet was used.'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    keepIndices = (
        parseRangeString(args.keepIndices, convertToZeroBased=True)
        if args.keepIndices else None)

    if args.keepIndicesFile:
        keepIndices = keepIndices or set()
        with open(args.keepIndicesFile) as fp:
            for lineNumber, line in enumerate(fp):
                try:
                    keepIndices.update(
                        parseRangeString(line, convertToZeroBased=True))
                except ValueError as e:
                    raise ValueError(
                        'Keep indices file %r line %d could not be parsed: %s'
                        % (args.keepIndicesFile, lineNumber, e))

    removeIndices = (
        parseRangeString(args.removeIndices, convertToZeroBased=True)
        if args.removeIndices else None)

    if args.removeIndicesFile:
        removeIndices = removeIndices or set()
        with open(args.removeIndicesFile) as fp:
            for lineNumber, line in enumerate(fp):
                try:
                    removeIndices.update(
                        parseRangeString(line, convertToZeroBased=True))
                except ValueError as e:
                    raise ValueError(
                        'Remove indices file %r line %d parse error: %s'
                        % (args.removeIndicesFile, lineNumber, e))

    reads.filter(
        minLength=args.minLength, maxLength=args.maxLength,
        removeGaps=args.removeGaps,
        whitelist=set(args.whitelist) if args.whitelist else None,
        blacklist=set(args.blacklist) if args.blacklist else None,
        whitelistFile=args.whitelistFile, blacklistFile=args.blacklistFile,
        titleRegex=args.titleRegex,
        negativeTitleRegex=args.negativeTitleRegex,
        truncateTitlesAfter=args.truncateTitlesAfter,
        indices=set(args.indices) if args.indices else None,
        head=args.head, removeDuplicates=args.removeDuplicates,
        removeDuplicatesById=args.removeDuplicatesById,
        removeDescriptions=args.removeDescriptions,
        randomSubset=args.randomSubset, trueLength=args.trueLength,
        sampleFraction=args.sampleFraction,
        sequenceNumbersFile=args.sequenceNumbersFile,
        keepIndices=keepIndices, removeIndices=removeIndices)

    saveAs = (
        args.saveAs or
        (args.fasta and 'fasta') or
        (args.fastq and 'fastq') or
        (args.fasta_ss and 'fasta-ss'))

    # Check for incompatible read/write formats. We can't write FASTQ
    # unless we have FASTQ on input (else we won't have quality information),
    # and we can't write PDB FASTA with secondary structure information
    # unless we have that on input.
    if saveAs == 'fastq' and not args.fastq:
        raise ValueError(
            'You have specified --saveAs fastq without using --fastq '
            'to indicate that the input is FASTQ. Please be explicit.')
    elif saveAs == 'fasta-ss' and not args.fasta_ss:
        raise ValueError(
            'You have specified --saveAs fasta-ss without using --fasta-ss '
            'to indicate that the input is PDB FASTA. Please be explicit.')

    write = sys.stdout.write
    kept = 0
    for read in reads:
        kept += 1
        write(read.toString(format_=saveAs))

    total = reads.unfilteredLength()

    if not args.quiet:
        print('Read %d sequence%s, kept %d (%.2f%%).' %
              (total, '' if total == 1 else 's', kept,
               0.0 if total == 0 else kept / total * 100.0), file=sys.stderr)

    if args.checkResultCount is not None:
        if kept != args.checkResultCount:
            if not args.quiet:
                print('Did not write the expected %d sequence%s (wrote %d).' %
                      (args.checkResultCount,
                       '' if args.checkResultCount == 1 else 's', kept),
                      file=sys.stderr)
            sys.exit(1)

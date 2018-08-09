#!/usr/bin/env python

from __future__ import print_function, division

import sys

from dark.filter import (
    addFASTAFilteringCommandLineOptions, parseFASTAFilteringCommandLineOptions,
    addFASTAEditingCommandLineOptions, parseFASTAEditingCommandLineOptions)
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


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
        '--checkResultCount', type=int,
        help=('The number of reads expected in the output. If this number is '
              'not seen, the script exits with status 1 and an error '
              'message is printed unless --quiet was used.'))

    addFASTACommandLineOptions(parser)
    addFASTAFilteringCommandLineOptions(parser)
    addFASTAEditingCommandLineOptions(parser)

    args = parser.parse_args()

    reads = parseFASTAEditingCommandLineOptions(
        args, parseFASTAFilteringCommandLineOptions(
            args, parseFASTACommandLineOptions(args)))

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

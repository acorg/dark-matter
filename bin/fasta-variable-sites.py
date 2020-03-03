#!/usr/bin/env python

import sys

from dark.dna import FloatBaseCounts
from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions, Reads)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin, write reads with only the '
                     'variable sites to stdout.'))

    parser.add_argument(
        '--printSites', action='store_true', default=False,
        help='Print the variable site locations to standad error.')

    parser.add_argument(
        '--gapChars', default='-',
        help='The characters that should be considered as gaps.')

    parser.add_argument(
        '--homogeneous', default=1.0, type=float,
        help=('If the most-common nucleotide frequency at a site is at least '
              'this value the site will be considered homogeneous.'))

    parser.add_argument(
        '--confirmed', action='store_true', default=False,
        help=('Only keep sites where there is confirmed variation (i.e., '
              'ambiguous sites that are compatible with there being no '
              'variation are not included).'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = list(parseFASTACommandLineOptions(args))

    if not reads:
        sys.exit(0)

    if len(set(len(read.sequence) for read in reads)) > 1:
        print('Input sequences are not all the same length!', file=sys.stderr)
        sys.exit(1)

    gaps = set(args.gapChars)
    confirmed = args.confirmed
    level = args.homogeneous

    variableSites = set()

    for site in range(len(reads[0].sequence)):
        codes = set(read.sequence[site] for read in reads) - gaps

        if codes:
            counts = FloatBaseCounts(codes)
            if counts.variable(confirmed) and not counts.homogeneous(level):
                # if counts.variable(confirmed):
                variableSites.add(site)
                if args.printSites:
                    print('%d: %s' % (site + 1, counts), file=sys.stderr)
        else:
            # All sites had a gap. This probably shouldn't happen.
            print('WARNING: Site %d had gaps in all reads!' % site + 1,
                  file=sys.stderr)

    saveAs = (args.fasta and 'fasta') or (args.fastq and 'fastq')

    Reads(reads).filter(keepSites=variableSites).save(
        sys.stdout, format_=saveAs)

    if args.printSites:
        n = len(variableSites)
        print('%d site%s %svariable (threshold for homogeneity: %.3f).' %
              (n,
               ' was' if n == 1 else 's were',
               'confirmed ' if confirmed else '',
               level), file=sys.stderr)

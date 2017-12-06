#!/usr/bin/env python

import sys
import argparse

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions, Reads)

parser = argparse.ArgumentParser(
    description=('Write sorted FASTA/Q to stdout. Sorting is by sequence id, '
                 'then by sequence, then by quality (if FASTQ)'))

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)

Reads(sorted(reads)).save(sys.stdout, 'fastq' if args.fastq else 'fasta')

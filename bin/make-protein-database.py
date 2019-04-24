#!/usr/bin/env python

from __future__ import print_function, division

import sys
import os
import warnings
from time import time
from itertools import chain
import argparse

from dark.proteins import SqliteIndex


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Create an sqlite3 database from AA sequences. '
                 'The protein sequences for the database are printed to '
                 'standard output.'))

parser.add_argument(
    '--databaseFile', required=True,
    help=('The output file. This file must not exist (use --force to '
          'overwrite).'))

parser.add_argument(
    '--force', default=False, action='store_true',
    help='If True and the output file already exists, overwrite it.')

parser.add_argument(
    '--quiet', default=False, action='store_true',
    help='Do not print indexing progress.')

parser.add_argument(
    '--noWarnings', default=False, action='store_true',
    help='Do not print warnings about unparseable GenBank records.')

parser.add_argument(
    '--gb', metavar='GenBank-file', nargs='+', action='append',
    required=True,
    help=('The GenBank file(s) to make the database from. These may be '
          'uncompressed, or compressed with bgzip (from samtools), with '
          'a .gz suffix.'))

parser.add_argument(
    '--accessionToNameFile', required=True,
    help=('The name of a file containing accession numbers and the '
          'description (i.e., full FASTA id line, sans leading >) of the '
          'corresponding nucleotide sequence.'))

parser.add_argument(
    '--nucleotideAccessionToTaxidFile', default='nucl_gb.accession2taxid.gz',
    help=('The nucleotide accession number to taxid file to read.  Must be '
          'gzipped (see ../doc/protein-database.md for detail).'))

args = parser.parse_args()

if os.path.exists(args.databaseFile):
    if args.force:
        os.unlink(args.databaseFile)
    else:
        print("Output file '%s' already exists. Use --force to overwrite."
              % args.databaseFile, file=sys.stderr)
        sys.exit(1)


def main(args):
    """
    Build the protein database.

    @param args: The namespace of command-line arguments returned by
        argparse.parse_args()
    """

    # Flatten the lists of lists that we get from using both nargs='+' and
    # action='append'. We use both because it allows people to use (e.g.)  --gb
    # on the command line either via "--gb file1 --gb file2" or "--gb file1
    # file2", or a combination of these. That way it's not necessary to
    # remember which way you're supposed to use it and you also can't be hit by
    # the subtle problem encountered in
    # https://github.com/acorg/dark-matter/issues/453
    gbFiles = list(chain.from_iterable(args.gb))

    verbose = not args.quiet

    # Read in the accession number to nucleotide sequence name file (this could
    # also have been obtained by just reading the full FASTA file of the
    # nucleotide sequences whose proteins are looked up).
    accessionToName = {}
    for line in open(args.accessionToNameFile):
        name, description = line[:-1].split('|', 1)
        assert name not in accessionToName
        accessionToName[name] = description

    if verbose:
        overallStart = time()

    with SqliteIndex(args.databaseFile) as db:
        for filename in gbFiles:
            if verbose:
                print("Indexing '%s' ... " % filename, end='', file=sys.stderr)
                start = time()

            genomeCount, proteinCount = db.addFile(filename, accessionToName)

            if verbose:
                elapsed = time() - start
                print('indexed %d sequence%s containing %d protein%s '
                      'in %.2f seconds.' %
                      (genomeCount, '' if genomeCount == 1 else 's',
                       proteinCount, '' if proteinCount == 1 else 's',
                       elapsed),
                      file=sys.stderr)

        db.updateGenomeTaxids(args.nucleotideAccessionToTaxidFile,
                              progressFp=(sys.stderr if verbose else None))
    if verbose:
        elapsed = time() - overallStart
        print('Overall indexing time: %.2f seconds (%.2f mins).' %
              (elapsed, elapsed / 60), file=sys.stderr)


if args.noWarnings:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        main(args)
else:
    main(args)

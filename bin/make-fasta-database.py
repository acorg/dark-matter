#!/usr/bin/env python

from __future__ import print_function

import sys
import os
from time import time
from itertools import chain

from dark.fasta import SqliteIndex


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Create an sqlite3 database from FASTA sequences.')

    parser.add_argument(
        '--out', required=True,
        help=('The output file. This file must not exist (use --force to '
              'overwrite).'))

    parser.add_argument(
        '--force', default=False, action='store_true',
        help='If True and the output file already exists, overwrite it.')

    parser.add_argument(
        '--quiet', default=False, action='store_true',
        help='If True do not print indexing progress.')

    parser.add_argument(
        '--fasta', metavar='FASTA-file', nargs='+', action='append',
        required=True,
        help=('the FASTA file(s) to make the database from. These may be '
              'uncompressed, or compressed with bgzip (from samtools), with '
              'a .gz suffix.'))

    args = parser.parse_args()

    if os.path.exists(args.out):
        if args.force:
            os.unlink(args.out)
        else:
            print("Output file '%s' already exists. Use --force to overwrite."
                  % args.out, file=sys.stderr)
            sys.exit(1)

    index = SqliteIndex(args.out)

    # Flatten the lists of lists that we get from using both nargs='+' and
    # action='append'. We use both because it allows people to use (e.g.)
    # --fasta on the command line either via "--fasta file1 --fasta file2"
    # or "--fasta file1 file2", or a combination of these. That way it's
    # not necessary to remember which way you're supposed to use it and you
    # also can't be hit by the subtle problem encountered in
    # https://github.com/acorg/dark-matter/issues/453
    fastaFiles = list(chain.from_iterable(args.fasta))

    verbose = not args.quiet

    for filename in fastaFiles:
        if verbose:
            print("Indexing '%s' ... " % filename, end='', file=sys.stderr)
            start = time()

        count = index.addFile(filename)

        if verbose:
            elapsed = time() - start
            print('indexed %d sequence%s in %.2f seconds.' %
                  (count, '' if count == 1 else 's', elapsed), file=sys.stderr)

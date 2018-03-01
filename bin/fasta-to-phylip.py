#!/usr/bin/env python

from __future__ import print_function

import sys
import tempfile

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA sequences (of equal length) on stdin '
                     'write Phylip to stdout.'))

    parser.add_argument(
        '--addIdSpace', action='store_true', default=False,
        help=('If True, print an extra space after each sequence id (the '
              'space is expected by some programs that process Phylip files, '
              'e.g., baseml).'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    idSpace = ' ' if args.addIdSpace else ''
    reads = parseFASTACommandLineOptions(args)

    count = 0
    length = None
    with tempfile.TemporaryFile() as fp:
        for read in reads:
            count += 1
            if length is None:
                length = len(read)
            else:
                if len(read) != length:
                    raise ValueError(
                        'FASTA sequence %r was not of the '
                        'expected length (%d)' % (read.id, length))

            name = read.id.split()[0]
            fp.write(('%s%s %s\n' %
                      (name, idSpace, read.sequence)).encode('utf-8'))

        if count:
            # Print the phylip header with the number of sequences and their
            # common length.
            print('%d %d' % (count, length))
            fp.seek(0)
            write = sys.stdout.write
            while True:
                chunk = fp.read(1024)
                if chunk:
                    write(chunk.decode('utf-8'))
                else:
                    break

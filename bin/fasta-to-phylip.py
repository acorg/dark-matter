#!/usr/bin/env python

from __future__ import print_function

import sys
import tempfile

from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.fastq import FastqReads


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin write the sequence ids and lengths '
                     'to stdout.'))

    parser.add_argument(
        '--readClass', default='fasta', choices=('fasta', 'fastq', 'fasta-ss'),
        help='If specified, give the input FASTA type')

    args = parser.parse_args()

    if args.readClass == 'fastq':
        # TODO: FastqReads should take a checkAlphabet argument, in the way
        # that FastaReads does.
        reads = FastqReads(sys.stdin)
    elif args.readClass == 'fasta':
        reads = FastaReads(sys.stdin, checkAlphabet=0)
    else:
        # args.readClass must be fasta-ss due to the 'choices' argument
        # passed to parser.add_argument value above.
        assert args.readClass == 'fasta-ss'
        reads = SSFastaReads(sys.stdin, checkAlphabet=0)

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
            fp.write(('%s %s\n' % (name, read.sequence)).encode('utf-8'))

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

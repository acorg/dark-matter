#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from dark import subtract

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate a fastafile, containing a specified '
        'set of reads',
        epilog='Given a list of readIds and a fastaFile, generate '
        'a new fastaFile based on these reads'
    )

    parser.add_argument(
        '--fasta', metavar='FASTA-files', type=str,
        help='a string of names of fastafiles which should be '
        'subtracted, deliminited by commas.')

    parser.add_argument(
        '--invert', default=False, action='store_true',
        help='If True, all reads in readIds will be returned, else, '
        'all reads not in readIds will be returned.')

    args = parser.parse_args()

    files = args.fasta
    fileList = files.split(', ')
    reads = subtract.fastaSubtract(fileList, invert=args.invert)

    SeqIO.write(reads, sys.stdout, 'fasta')
    print >>sys.stderr, 'Found %d sequences.' % len(reads)

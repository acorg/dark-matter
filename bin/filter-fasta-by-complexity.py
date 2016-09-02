#!/usr/bin/env python

"""
Read DNA FASTA from stdin and print FASTA to stdout, excluding sequences
with too many bases in low complexity regions.

You must have dustmasker in your shell's PATH for this code to succeed.
dustmasker is part of the NCBI BLAST+ suite. You can get it from
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
"""

from __future__ import print_function

import sys
from os import close, unlink
from time import time
import argparse
from subprocess import call
from tempfile import mkstemp

from dark.fasta import FastaReads


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter FASTA for low complexity sequences',
        epilog='Read DNA FASTA from stdin and print FASTA to stdout, '
        'excluding sequences of low complexity.')

    parser.add_argument(
        'fasta', help='The FASTA file to read.')

    parser.add_argument(
        '--maxLowComplexity', metavar='FRAC', type=float, default=0.1,
        help='Sequences with a fraction of low-complexity bases up to this '
        'level will be filtered out.')

    parser.add_argument(
        '--rejectFile', metavar='FILE', default=None,
        help='The name of a file to save rejected sequences to.')

    args = parser.parse_args()
    fd, outFilename = mkstemp()
    close(fd)

    if args.rejectFile is not None:
        rejectFp = open(args.rejectFile, 'w')
        reject = rejectFp.write
    else:
        reject = lambda x: None

    cmd = 'dustmasker -in "%s" -out "%s" -outfmt fasta' % (
        args.fasta, outFilename)
    print('Running dustmasker. Hang in there... ', end=' ', file=sys.stderr)
    start = time()

    try:
        status = call(cmd, shell=True)
        print('done.\ndustmasker completed in %.2f seconds.' % (
            time() - start), file=sys.stderr)
        if status < 0:
            # Child terminated by signal.
            unlink(outFilename)
            sys.exit(-status)
        elif status > 0:
            print('WARNING: Child exited with status', status, file=sys.stderr)
            print('Output left in', outFilename, file=sys.stderr)
            sys.exit(status)
    except OSError as e:
        print('Execution failed:', e, file=sys.stderr)
        unlink(outFilename)
        sys.exit(1)

    save = sys.stdout.write
    reads = FastaReads(outFilename)
    maxLowComplexity = args.maxLowComplexity
    readCount = rejectCount = 0

    for read in reads:
        readCount += 1
        fasta = read.toString('fasta')
        if read.lowComplexityFraction() < maxLowComplexity:
            save(fasta)
        else:
            reject(fasta)
            rejectCount += 1

    unlink(outFilename)

    if args.rejectFile is not None:
        rejectFp.close()

    print((
        '%d sequences read, %d (%.2f%%) saved, %d (%.2f%%) rejected.' % (
            readCount, readCount - rejectCount,
            (readCount - rejectCount) / float(readCount) * 100.0,
            rejectCount,
            rejectCount / float(readCount) * 100.0)), file=sys.stderr)

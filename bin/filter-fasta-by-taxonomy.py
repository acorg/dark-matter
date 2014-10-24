#!/usr/bin/env python

"""
Read DNA FASTA from stdin and print FASTA to stdout, only including sequences
that match the requested taxonomy regular expression at any level.

This script will only produce output if the FASTA it is given has description
lines like

  >gi|9630656|ref|NC_001925.1| Bunyamwera virus L segment, complete sequence

with a GI number in the second '|'-separated field. This is the format found
in the NCBI reference sequence files. E.g. at
ftp://ftp.ncbi.nih.gov/refseq/release/viral

You must have the NCBI taxonomy database installed locally for this code to
succeed.  See the file doc/taxonomy.md for instructions on how to set that up.
If you see sequences unexpectedly rejected because they have no associated
taxonomy, make sure you have the latest taxonomy files loaded into MySQL.
"""

import sys
import argparse
import re

from dark.fasta import FastaReads
from dark.taxonomy import LineageFetcher


def writeDetails(accept, readId, taxonomy, fp):
    """
    Write read and taxonomy details.

    @param accept: A C{bool} indicating whether the read was accepted,
        according to its taxonomy.
    @param readId: The C{str} id of the read.
    @taxonomy: A C{list} of taxonomy C{str} levels.
    @fp: An open file pointer to write to.
    """
    fp.write('%s: %s\n        %s\n\n' % (
        'ACCEPT' if accept else 'REJECT', readId,
        ' | '.join(taxonomy) if taxonomy else 'No taxonomy found.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter FASTA based on taxonomy',
        epilog='Read DNA FASTA from stdin and print FASTA to stdout, only '
        'including sequences that match the requested taxonomy regular '
        'expression at any level.')

    parser.add_argument(
        '--taxonomy', type=str, required=True,
        help='The regex to match the taxonomy on. Case is ignored.')

    parser.add_argument(
        '--rejects', metavar='FILE', type=str, default=None,
        help='The name of a file to save rejected sequences to.')

    parser.add_argument(
        '--details', metavar='FILE', type=str, default=None,
        help='The name of a file to save taxonomy details to')

    args = parser.parse_args()

    if args.rejects is not None:
        rejectFp = open(args.rejects, 'w')
        reject = rejectFp.write
    else:
        rejectFp = None
        reject = lambda x: None

    if args.details is not None:
        detailsFp = open(args.details, 'w')
        details = lambda accept, readId, taxonomy: (
            writeDetails(accept, readId, taxonomy, detailsFp))
    else:
        detailsFp = None
        details = lambda accept, readId, taxonomy: None

    lineageFetcher = LineageFetcher()
    reads = FastaReads(sys.stdin)
    regexp = re.compile(args.taxonomy, re.I)
    save = sys.stdout.write
    readCount = rejectCount = 0

    for read in reads:
        readCount += 1
        fasta = read.toString('fasta')
        taxonomy = lineageFetcher.lineage(read.id)
        accept = False
        for level in taxonomy:
            if regexp.match(level):
                save(fasta)
                accept = True
                break
        else:
            rejectCount += 1
            reject(fasta)

        details(accept, read.id, taxonomy)

    if rejectFp:
        rejectFp.close()

    if detailsFp:
        detailsFp.close()

    lineageFetcher.close()

    print >>sys.stderr, (
        '%d sequences read, %d (%.2f%%) saved, %d (%.2f%%) rejected.' % (
            readCount, readCount - rejectCount,
            (readCount - rejectCount) / float(readCount) * 100.0,
            rejectCount,
            rejectCount / float(readCount) * 100.0))

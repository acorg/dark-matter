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
from typing import Optional, TextIO

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
    fp.write(
        "%s %s\n       %s\n\n"
        % (
            "MATCH:" if accept else "MISS: ",
            readId,
            " | ".join(taxonomy) if taxonomy else "No taxonomy found.",
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter FASTA based on taxonomy",
        epilog="Read DNA FASTA from stdin and print FASTA to stdout, only "
        "including sequences that match the requested taxonomy regular "
        "expression at any level.",
    )

    parser.add_argument(
        "--taxonomy",
        required=True,
        help="The regex to match the taxonomy on. Case is ignored.",
    )

    parser.add_argument(
        "--invert",
        action="store_true",
        default=False,
        help="If True, only write sequences whose taxonomy does not match.",
    )

    parser.add_argument(
        "--detailsFile",
        metavar="FILE",
        default=None,
        help="The name of a file to save taxonomy details to",
    )

    args = parser.parse_args()

    try:
        regexp = re.compile(args.taxonomy, re.I)
    except re.error as e:
        print(
            "Could not compile %r to a regular expression:" % args.taxonomy,
            e,
            file=sys.stderr,
        )
        sys.exit(1)

    detailsFp: Optional[TextIO]

    if args.detailsFile is not None:
        detailsFp = open(args.detailsFile, "w")

        def details(accept, readId, taxonomy):
            return writeDetails(accept, readId, taxonomy, detailsFp)

    else:
        detailsFp = None

        def details(accept, readId, taxonomy):
            return None

    lineageFetcher = LineageFetcher()
    reads = FastaReads(sys.stdin)
    save = sys.stdout.write
    readCount = saveCount = noTaxonomyCount = 0

    for read in reads:
        readCount += 1
        fasta = read.toString("fasta")
        taxonomy = lineageFetcher.lineage(read.id)
        if taxonomy:
            for taxonomyId, scientificName in taxonomy:
                if regexp.match(scientificName):
                    details(True, read.id, taxonomy)
                    if not args.invert:
                        saveCount += 1
                        save(fasta)
                    break
            else:
                details(False, read.id, taxonomy)
                if args.invert:
                    saveCount += 1
                    save(fasta)
        else:
            noTaxonomyCount += 1

    if detailsFp:
        detailsFp.close()

    lineageFetcher.close()

    rejectCount = readCount - saveCount - noTaxonomyCount
    print(
        "%d sequences read, %d (%.2f%%) saved, %d (%.2f%%) rejected, "
        "%d (%.2f%%) no taxomony found."
        % (
            readCount,
            saveCount,
            saveCount / float(readCount) * 100.0,
            rejectCount,
            (rejectCount) / float(readCount) * 100.0,
            noTaxonomyCount,
            noTaxonomyCount / float(readCount) * 100.0,
        ),
        file=sys.stderr,
    )

#!/usr/bin/env python

"""
Verify that a FASTA file and our JSON BLAST file that was produced from it
do not violate some basic compatibility checks.

Exits non-zero with an error message on stderr if there's a problem, else
exits silently with zero status.
"""

import sys

from dark.blast.alignments import BlastReadsAlignments
from dark.fasta import FastaReads


def check(fastaFile, jsonFiles):
    """
    Check for simple consistency between the FASTA file and the JSON files.

    Note that some checking is already performed by the BlastReadsAlignments
    class. That includes checking the number of reads matches the number of
    BLAST records and that read ids and BLAST record read ids match.

    @param jsonFiles: A C{list} of names of our BLAST JSON. These may
        may be compressed (as bz2).
    @param fastaFile: The C{str} name of a FASTA-containing file.
    """
    reads = FastaReads(fastaFile)
    readsAlignments = BlastReadsAlignments(reads, jsonFiles)
    for index, readAlignments in enumerate(readsAlignments):
        # Check that all the alignments in the BLAST JSON do not have query
        # sequences or query offsets that are greater than the length of
        # the sequence given in the FASTA file.
        fastaLen = len(readAlignments.read)
        for readAlignment in readAlignments:
            for hsp in readAlignment.hsps:
                # The FASTA sequence should be at least as long as the
                # query in the JSON BLAST record (minus any gaps).
                assert fastaLen >= len(hsp.query) - hsp.query.count("-"), (
                    "record %d: FASTA len %d < HSP query len %d.\n"
                    "FASTA: %s\nQuery match: %s"
                    % (
                        index,
                        fastaLen,
                        len(hsp.query),
                        readAlignments.read.sequence,
                        hsp.query,
                    )
                )
                # The FASTA sequence length should be larger than either of
                # the query offsets mentioned in the JSON BLAST
                # record. That's because readStart and readEnd are offsets
                # into the read - so they can't be bigger than the read
                # length.
                #
                # TODO: These asserts should be more informative when they
                # fail.
                assert fastaLen >= hsp.readEnd >= hsp.readStart, (
                    "record %d: FASTA len %d not greater than both read "
                    "offsets (%d - %d), or read offsets are non-increasing. "
                    "FASTA: %s\nQuery match: %s"
                    % (
                        index,
                        fastaLen,
                        hsp.readStart,
                        hsp.readEnd,
                        readAlignments.read.sequence,
                        hsp.query,
                    )
                )


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        check(sys.argv[1], sys.argv[2:])
    else:
        print(
            "Usage: %s file.fasta file1.json file2.json..." % sys.argv[0],
            file=sys.stderr,
        )
        sys.exit(1)

#!/usr/bin/env python

import sys

from prseq import FastaReader

from dark.reads import Reads
from dark.sequence import findPrimerBidiLimits


def trimPrimers(primer: str, verbose: bool) -> None:
    """
    @param primer: A string sequence.
    @param verbose: A C{bool}, if C{True} output additional information about
        how often and where primers were found.
    """
    reads = Reads()
    absentCount = forwardCount = reverseCount = count = 0
    for seqRecord in FastaReader():
        count += 1
        start, end = findPrimerBidiLimits(primer, seqRecord.sequence)
        if start == 0:
            if end == len(seqRecord):
                absentCount += 1
            else:
                reverseCount += 1
        else:
            forwardCount += 1
            if end != len(seqRecord):
                reverseCount += 1
        reads.add(seqRecord[start:end])

    if verbose:
        print(
            ("Read %d sequences. Found forward: %d, " "Found reversed: %d, Absent: %d")
            % (count, forwardCount, reverseCount, absentCount),
            file=sys.stderr,
        )

    reads.save(sys.stdout)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Given FASTA on stdin, look for a primer sequence "
            "and write trimmed FASTA (after the primer) to stdout."
        )
    )

    parser.add_argument("primer", help="the primer sequence")
    parser.add_argument(
        "--verbose",
        type=bool,
        help="If True, print information on found primers",
    )

    args = parser.parse_args()

    trimPrimers(args.primer.upper(), args.verbose)

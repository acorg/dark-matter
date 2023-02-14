#!/usr/bin/env python

"""
Aligns each FASTA sequence in one file to each sequence in another.
"""

import argparse

from dark.fasta import FastaReads
from dark.local_align import LocalAlignment


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Given two FASTA files, perform a local Smith Waterman alignment "
        "on all the sequences in one file against all those in the other."
    ),
)

parser.add_argument("--fastaFile1", required=True, help="The first FASTA file.")

parser.add_argument("--fastaFile2", required=True, help="The second FASTA file.")

parser.add_argument("--matchScore", type=int, default=1, help="The match score.")

parser.add_argument("--mismatchScore", type=int, default=-1, help="The mismatch score.")

parser.add_argument("--gapOpenScore", type=int, default=-1, help="The gap open score.")

parser.add_argument(
    "--gapExtendScore", type=int, default=-1, help="The gap extend score."
)

parser.add_argument(
    "--gapExtendDecay", type=float, default=0.0, help="The gap extend decay."
)

args = parser.parse_args()

for seq1 in FastaReads(args.fastaFile1):
    for seq2 in FastaReads(args.fastaFile2):
        alignment = LocalAlignment(
            seq1,
            seq2,
            match=args.matchScore,
            mismatch=args.mismatchScore,
            gap=args.gapOpenScore,
            gapExtend=args.gapExtendScore,
            gapExtendDecay=args.gapExtendDecay,
        )

        print(alignment.createAlignment(resultFormat=str))

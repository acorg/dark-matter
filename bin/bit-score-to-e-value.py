#!/usr/bin/env python

"""
Convert a bit score to an e-value.

Example usage:

$ bit-score-to-evalue.py --dbSize 168142520 --dbSequenceCount 5660 \
  --queryLength 111 --lengthAdjustment 26 --bitScore 37.3537
0.0813077725194
"""

import argparse

from dark.blast.score import bitScoreToEValue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a bit score to an e-value.")

    parser.add_argument(
        "--bitScore", type=float, required=True, help="The bit score to convert"
    )

    parser.add_argument(
        "--dbSize",
        type=int,
        required=True,
        help="The total number of bases in the sequence database.",
    )

    parser.add_argument(
        "--dbSequenceCount",
        type=int,
        required=True,
        help="The number of sequences in the database.",
    )

    parser.add_argument(
        "--queryLength",
        type=int,
        required=True,
        help="The length of the query sequence.",
    )

    parser.add_argument(
        "--lengthAdjustment", type=int, required=True, help="The length adjustment."
    )

    args = parser.parse_args()
    print(
        bitScoreToEValue(
            args.bitScore,
            args.dbSize,
            args.dbSequenceCount,
            args.queryLength,
            args.lengthAdjustment,
        )
    )

#!/usr/bin/env python

"""
Convert an e-value to a bit score.

Example usage:

$ e-value-to-bit-score.py --dbSize 168142520 --dbSequenceCount 5660 \
  --queryLength 111 --lengthAdjustment 26 --eValue 0.0813077725194
37.3537
"""

from __future__ import print_function

import argparse

from dark.blast.score import eValueToBitScore


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert a bit score to an e-value.')

    parser.add_argument(
        '--eValue', type=float, required=True, help='The e-value to convert.')

    parser.add_argument(
        '--dbSize', type=int, required=True,
        help='The total number of bases in the sequence database.')

    parser.add_argument(
        '--dbSequenceCount', type=int, required=True,
        help='The number of sequences in the database.')

    parser.add_argument(
        '--queryLength', type=int, required=True,
        help='The length of the query sequence.')

    parser.add_argument(
        '--lengthAdjustment', type=int, required=True,
        help='The length adjustment.')

    args = parser.parse_args()
    print(eValueToBitScore(args.eValue, args.dbSize, args.dbSequenceCount,
                           args.queryLength, args.lengthAdjustment))

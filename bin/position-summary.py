#!/usr/bin/env python

import argparse

from dark.fasta import FastaReads


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "For a fasta file with sequences, summarize what is "
            "happening at a specific position."
        ),
    )

    parser.add_argument(
        "--fastaFile", required=True, help="The name of the FASTA file to read."
    )

    parser.add_argument(
        "--position",
        required=True,
        type=int,
        help="The (one-based) position to examine.",
    )

    args = parser.parse_args()
    reads = FastaReads(args.fastaFile)
    result = reads.summarizePosition(args.position - 1)
    nReads = result["readCount"]

    print(
        "%d of %d sequences were excluded due to length."
        % (result["excludedCount"], nReads)
    )

    denominator = (nReads - result["excludedCount"]) / 100.0
    for base, count in result["countAtPosition"].items():
        print("%s: Total: %s (%.2f%%)" % (base, count, count / denominator))

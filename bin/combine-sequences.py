#!/usr/bin/env python

import argparse
import multiprocessing
import sys

from dark.aligners import mafft
from dark.fasta import FastaReads
from dark.reads import DNARead

MAFFT_DEFAULT_ARGS = "--globalpair --maxiterate 1000 --preservecase"
MAFFT_ALGORITHMS_URL = (
    "https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html"
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Combine multiple sequences into a single sequence.",
    )

    parser.add_argument(
        "--fastaFiles",
        required=True,
        nargs="+",
        help=("The name of the FASTA files containing the sequences to be combined."),
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help=(
            "The number of threads to use when running the aligner (if "
            "--align is used and the alignment algorithm can make use of "
            "multiple threads)."
        ),
    )

    parser.add_argument(
        "--alignerOptions",
        default=MAFFT_DEFAULT_ARGS,
        help=(
            "Optional arguments to pass to the alignment algorithm. The "
            "default options are %r. See %s for some possible option "
            "combinations." % (MAFFT_DEFAULT_ARGS, MAFFT_ALGORITHMS_URL)
        ),
    )

    args = parser.parse_args()

    reads = []
    for fastaFile in args.fastaFiles:
        reads.extend(FastaReads(fastaFile))

    if len(reads) == 1:
        print("Cannot combine just one read. Exiting.", file=sys.stderr)
        sys.exit(1)

    alignedReads = mafft(reads, options=args.alignerOptions, threads=args.threads)

    combinedSequence = alignedReads.combineReads()

    print(DNARead("Combined", combinedSequence).toString(format_="fasta"))

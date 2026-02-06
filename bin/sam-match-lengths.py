import argparse
import csv
import sys
from contextlib import redirect_stdout
from pathlib import Path
from typing import Sequence

import numpy as np
import pysam

from dark.utils import pct


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Print information about the lengths of query matches in a SAM/BAM file. "
            "If --tsv is given, print (to standard output) tab-separated fields "
            "containing: the query id, query lengths, and matched region length (i.e., "
            "excluding soft-clipped bases). The number of soft-clipped bases on the "
            "left and right are also given in the last two TSV fields. Queries that "
            "did not align produce no output. If --stats is given, summary statistics "
            "are printed to standard output (or to standard error if --tsv is also "
            "given)."
        ),
    )

    parser.add_argument(
        "samfile",
        type=Path,
        help="The SAM/BAM file to read.",
    )

    parser.add_argument(
        "--tsv",
        action="store_true",
        help="Print TSV to standard output.",
    )

    parser.add_argument(
        "--stats",
        action="store_true",
        help=(
            "Print summary statistics to standard output (or standard error if --tsv "
            "is given)."
        ),
    )

    return parser.parse_args()


def summarize(what: str, lengths: Sequence[int]) -> None:
    min_ = np.min(lengths)
    max_ = np.max(lengths)
    mean = np.mean(lengths)
    median = np.median(lengths)
    sd = np.std(lengths)

    print(f"{what}:")
    print(f"  N: {len(lengths):,}")
    print(f"  Min: {min_:,}")
    print(f"  Max: {max_:,}")
    print(f"  Mean: {mean:.3f}")
    print(f"  Median: {median:,}")
    print(f"  SD: {sd:.3f}")


def main() -> None:
    args = get_args()
    tsv = args.tsv
    stats = args.stats

    if not (tsv or stats):
        stats = True

    af = pysam.AlignmentFile(str(args.samfile))
    CSOFT_CLIP = pysam.CIGAR_OPS.CSOFT_CLIP

    if tsv:
        writerow = csv.writer(sys.stdout, delimiter="\t").writerow
        writerow(
            (
                "Query",
                "Length",
                "Soft-clipped length",
                "Soft-clipped left",
                "Soft-clipped right",
            )
        )

    if stats:
        lengths = []
        clipped_lengths = []
        clipped_lengths_left = []
        clipped_lengths_right = []
        read_count = mapped_count = 0

    for alignment in af.fetch():
        if stats:
            read_count += 1
        if alignment.is_mapped:
            if stats:
                mapped_count += 1
            ct = alignment.cigartuples
            soft_left = ct[0][1] if ct[0][0] == CSOFT_CLIP else 0
            soft_right = ct[-1][1] if len(ct) > 1 and ct[-1][0] == CSOFT_CLIP else 0
            length = alignment.query_length
            clipped_length = length - (soft_left + soft_right)

            if stats:
                lengths.append(length)
                clipped_lengths.append(clipped_length)
                clipped_lengths_left.append(soft_left)
                clipped_lengths_right.append(soft_right)

            if tsv:
                writerow(
                    (
                        alignment.query_name,
                        length,
                        clipped_length,
                        soft_left,
                        soft_right,
                    )
                )

    if stats:
        with redirect_stdout(sys.stderr if tsv else sys.stdout):
            print(f"Number of mapped reads: {pct(mapped_count, read_count)}.")
            total_clipped = sum(clipped_lengths_left) + sum(clipped_lengths_right)
            print(f"Total soft-clipped bases: {pct(total_clipped, sum(lengths))}.")
            summarize("Query lengths", lengths)
            summarize("Soft-trimmed lengths", clipped_lengths)
            summarize("Left-trimmed lengths", clipped_lengths_left)
            summarize("Right-trimmed lengths", clipped_lengths_right)


if __name__ == "__main__":
    main()

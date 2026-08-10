#!/usr/bin/env python
"""
Merge, sort, and index a set of SAM/BAM files using samtools.

Given a list of input SAM/BAM files, this script:
  1. Concatenates them into a single BAM file (samtools cat)
  2. Sorts the result by coordinate (samtools sort)
  3. Indexes the sorted BAM (samtools index) to produce a .bai file

Requires samtools to be installed and available on PATH.

Original version written by Claude on 2026-08-10.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate, sort, and index SAM/BAM files with samtools."
    )
    parser.add_argument(
        "--bam",
        required=True,
        type=Path,
        help="Output BAM filename.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        metavar="FILE",
        help="Input SAM/BAM files to merge.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )
    parser.add_argument(
        "-@",
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for sorting/indexing (default: 1).",
    )
    return parser.parse_args()


def check_samtools():
    if shutil.which("samtools") is None:
        sys.exit("Error: samtools not found on PATH. Please install samtools.")


def run(cmd):
    """Run a command, streaming output, and abort on failure."""
    cmd = [str(c) for c in cmd]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error: command {' '.join(cmd)!r} failed: {e}.")


def main():
    args = parse_args()
    check_samtools()

    for f in args.inputs:
        if not f.is_file():
            sys.exit(f"Error: input file {str(f)!r} not found.")

    if args.bam.exists() and not args.force:
        sys.exit(
            f"Error: output file '{args.bam}' already exists. "
            "Use --force/-f to overwrite."
        )

    out_dir = args.bam.resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=out_dir) as tmpdir:
        merged_path = Path(tmpdir) / "merged.bam"

        run(["samtools", "cat", "-o", merged_path, *args.inputs])

        run(
            [
                "samtools",
                "sort",
                "-@",
                args.threads,
                "-o",
                args.bam,
                merged_path,
            ]
        )

    run(["samtools", "index", "-@", args.threads, args.bam])


if __name__ == "__main__":
    main()

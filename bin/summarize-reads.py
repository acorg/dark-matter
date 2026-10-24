#!/usr/bin/env python

import argparse

from dark.summarize import summarizeReads

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Given a FASTA or FASTQ file, print summary statistics: read "
            "count, total/average/median/min/max length, and per-base "
            "counts."
        )
    )
    parser.add_argument("filename", help="The FASTA or FASTQ file to summarize.")
    args = parser.parse_args()

    filename = args.filename
    if filename.lower().endswith("q"):
        filetype = "fastq"
    else:
        # Assume anything not ending with a 'q' is fasta.
        filetype = "fasta"

    result = summarizeReads(filename, filetype)

    print("Number of reads:", result["read_number"])
    print("Total length: %s bases" % result["total_length"])
    print("The average read length: %s bases" % result["average_length"])
    print("Longest read:", result["max_length"])
    print("Shortest read:", result["min_length"])
    print("Median length:", result["median_length"])

    for base in sorted(result["base_counts"]):
        count = result["base_counts"][base]
        print(
            "%s: Total: %s; Average per read: %s"
            % (base, count, count / result["read_number"])
        )

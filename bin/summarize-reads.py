#!/usr/bin/env python

from dark.summarize import summarizeReads
import sys


if len(sys.argv) != 2:
    print("Usage: %s file.fasta / file.fastq" % sys.argv[0], file=sys.stderr)
    sys.exit(1)
else:
    filename = sys.argv[1]
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

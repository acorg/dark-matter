#!/usr/bin/env python

from dark.summarize import summarize_reads
import sys


if len(sys.argv) > 2:
    print >>sys.stderr, "Usage: %s file.fasta / file.fastq" % sys.argv[0]
    sys.exit(1)

else:
    filename = sys.argv[1]
    result = summarize_reads(filename)

    print "Number of reads:", result["read_number"]
    print "Total length: %s bases" % result["total_length"]
    print "The average read length: %s bases" % result["average_length"]
    print "Longest read:", result["max_length"]
    print "Shortest read:", result["min_length"]
    print "Median length:", result["median_length"]

    for base, count in result["base_counts"].items():
        print "%s: Total: %s; Average per read: %s" % (
            base, count, count / result["read_number"])

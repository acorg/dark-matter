#!/bin/sh

# Read SAM (see http://samtools.github.io/hts-specs/SAMv1.pdf) format from
# stdin and write FASTQ (https://en.wikipedia.org/wiki/FASTQ_format) to
# stdout.
#
# Note that the '+' line of each sequence in the FASTQ output also has the
# sequence id (which is present in the first line). Although this is
# supposedly optional (according to the Wikipedia link above), some tools
# that process FASTQ insist on it being present - so we leave it in.

egrep -v '^@' | awk '{printf "@%s\n%s\n+%s\n%s\n", $1, $10, $1, $11}'

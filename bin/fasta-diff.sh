#!/bin/bash

set -Eeuo pipefail

# This is a crude FASTA/Q diff. It takes advantage of the system sort and
# diff by converting FASTA/FASTQ records to a single TAB-separated line
# each, sorting them, then calling diff. The result is then converted back
# to FASTA/FASTQ (though with < and > indicators from diff preceeding
# sequence ids to show in which input file they were found). Line number
# output from the diff is discarded as it has no relation to the input
# files (the input sequences and qualities may be spread over many lines
# and the input sequences may be in any order, whereas the files passed to
# diff have ids, sequences, and quality strings on a single line and will
# be in sorted order so as to minimize diff output).

base=$(basename $0)

function oneLine()
{
    TMP=$(mktemp /tmp/${base}.XXXXXX) || {
        echo "$name: could not create temp file." >&2
        exit 2
    }

    case "$1" in
        *.fastq | *.FASTQ) bin/fasta-join.py --fastq < "$1" | sort > $TMP ;;
        *) bin/fasta-join.py < "$1" | sort > $TMP ;;
    esac

    echo $TMP
}

case $# in
    2)
        TMP1=$(oneLine "$1")
        TMP2=$(oneLine "$2")
        ;;
    *)
        echo "Usage: $base file1.fast[aq] file1.fast[aq]" >&2
        exit 1
        ;;
esac

# Convert the diff output lines back to FASTA/FASTQ.
#
# Use || true in the following to avoid exiting if diff finds a difference
# and exits non-zero, causing the script to exit due to the set -e.
{ diff "$TMP1" "$TMP2" | egrep -e '^[<>] ' | tr '\t' '\n'; } || true

rm "$TMP1" "$TMP2"

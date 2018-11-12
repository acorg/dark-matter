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

    catter=cat
    fastq=''

    case "$1" in
        *.fq | *.fastq | *.FASTQ) fastq='--fastq';;
        *.fq.gz | *.fastq.gz | *.FASTQ.gz) catter=gunzip; fastq='--fastq';;
        *.fq.bz2 | *.fastq.bz2 | *.FASTQ.bz2) catter=bzcat; fastq='--fastq';;
        *.gz | *.gz) catter=gunzip;;
        *.bz2 | *.bz2) catter=bzcat;;
    esac

    $catter < "$1" | fasta-join.py $fastq | sort > $TMP

    echo $TMP
}

diffArgs=

while [ $# -gt 0 ]
do
    case "$1" in
        -*) diffArgs="$diffArgs $1"; shift;;
        *) break;;
    esac
done

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
TMP3=$(mktemp /tmp/${base}.XXXXXX) || {
    echo "$name: could not create temp file." >&2
    exit 2
}

# Don't exit if the diff exits non-zero.
set +e
diff $diffArgs "$TMP1" "$TMP2" > "$TMP3"
status=$?
set -e

# Filter and emit the diff output.
egrep -e '^[<>] ' < "$TMP3" | tr '\t' '\n'

rm "$TMP1" "$TMP2" "$TMP3"

exit $status

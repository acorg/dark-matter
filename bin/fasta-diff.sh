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
#
# Thus this tool doesn't really help you find where differences are, it
# just shows you a concise summary of differences found. This can be very
# useful when there are no differences or when the number of differences is
# small or when the only differences are in the sequence ids (use
# --removeIds to remove them from the comparison). If the two input files
# are totally different, the output won't be of much use, except to
# indicate that there's a large difference.

base=$(basename $0)

function tempfile()
{
    tmp=$(mktemp /tmp/${base}.XXXXXX) || {
        echo "$name: could not create temp file." >&2
        exit 2
    }

    echo "$tmp"
}

function oneLine()
{
    # Echo a command line that will join a FASTA file (using fasta-join.py)
    # and sort the result.

    file="$1"
    tmp="$2"
    removeIds="$3"

    catter=stdin
    fastq=''

    case "$file" in
        *.fq | *.fastq | *.FASTQ) fastq='--fastq';;
        *.fq.gz | *.fastq.gz | *.FASTQ.gz) catter=gunzip; fastq='--fastq';;
        *.fq.bz2 | *.fastq.bz2 | *.FASTQ.bz2) catter=bzcat; fastq='--fastq';;
        *.gz | *.gz) catter=gunzip;;
        *.bz2 | *.bz2) catter=bzcat;;
    esac

    case $catter in
        stdin) echo "fasta-join.py $removeIds $fastq < \"$file\" | sort > \"$tmp\"";;
        *) echo "$catter < \"$file\" | fasta-join.py $removeIds $fastq | sort > \"$tmp\"";;
    esac
}

diffArgs=
removeIds=

while [ $# -gt 0 ]
do
    case "$1" in
        --removeIds)
            removeIds=--removeIds
            shift
            ;;
        -*)
            diffArgs="$diffArgs $1"
            shift
            ;;
        *)
            break
            ;;
    esac
done

case $# in
    2)
        tmp1=$(tempfile)
        tmp2=$(tempfile)
        command1=$(oneLine "$1" "$tmp1" "$removeIds")
        command2=$(oneLine "$2" "$tmp2" "$removeIds")

        # Run the setup commands in parallel if GNU parallel is installed.
        if [ $(type -t parallel) = 'filexx' ]
        then
            echo -e "$command1\n$command2" | parallel
        else
            eval $command1
            eval $command2
        fi
        ;;
    *)
        echo "Usage: $base [--removeIds] file1.fast[aq] file1.fast[aq]" >&2
        exit 1
        ;;
esac

# Run diff, keep its exit status, and convert its output lines back to
# FASTA/FASTQ.
tmp3=$(tempfile)

# Don't exit if the diff exits non-zero.
set +e
diff $diffArgs "$tmp1" "$tmp2" > "$tmp3"
status=$?
set -e

# Filter and emit the diff output. The tr will not do anything if
# --removeIds was used and the input was FASTA (not FASTQ), because in that
# case the output of fasta-join.py has no TABs in it. But let's not bother
# trying to detect that case.
#
# The egrep just finds the diff output lines that show differences (these
# are marked with < or >). This throws away the diff output lines numbers,
# but those were of no use anyway since we sorted the sequences in the
# FASTA/Q files.
egrep -e '^[<>] ' < "$tmp3" | tr '\t' '\n'

rm "$tmp1" "$tmp2" "$tmp3"

exit $status

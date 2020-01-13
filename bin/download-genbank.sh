#!/bin/bash

# This script downloads flat files (in GenBank format) of sequence info
# from NCBI. It is one of a set of utilities that we use to make a protein
# database.
#
# Usage is e.g.:
#
#     download-genbank.sh [-n] [id-dir] < accession-numbers-file
#
# where:
#
#     accession-numbers-file is a file containing GenBank accession numbers,
#         one per line. These should include the version number.
#
#     id-dir is the name of a directory in which to store files of sequence
#         ids, parameters to cURL, and GenBank files downloaded from
#         NCBI. If not given, defaults to the current directory.
#
#     -n option can be used to make the script just print out what it would
#         download, rather than actually doing the downloads. In this case,
#         a temporary id directory will be created and used and its name
#         will be printed on standard output.
#
# See ../doc/protein-database.md for information on the various tools
# needed to make a protein database and ../misc/Makefile-protein-database
# for a Makefile that can do the work of calling the scripts in order with
# a reasonable set of output file and directory names.

# Note that the number of [a-z] repeats in filenames in various places in
# the code below matches the 5 in the call to split.

set -Eeuo pipefail

URL=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
maxIdsPerRequest=200

# These are accession numbers that correspond to dead (suppressed or
# withdrawn) sequence records (see ../doc/protein-database.md) for more
# info. These were all in the RVDB clustered database in early 2019.
DEAD_REGEX='JX873962|KX912841|NC_021196|NC_036584'

dryRun=0

while [ $# -gt 0 ]
do
    case "$1" in
        -n)
            dryRun=1
            shift
            ;;
        *)
            break
            ;;
    esac
done

case $# in
    0) idDir=.;;
    1) idDir=$1;;
    *) echo "Usage: $(basename $0) [-n] [id-dir]" >&2; exit 1;;
esac

if [ $dryRun -eq 0 ]
then
    test -d $idDir || mkdir "$idDir"
else
    # Use a temporary id directory for a dry run.
    idDir=$(mktemp -d)
fi

# Read stdin for accession numbers and split them into files we'll use in
# entrez requests to NCBI.
egrep -v "$DEAD_REGEX" | split -a 5 -l $maxIdsPerRequest - $idDir/ids_
count=$(cat $idDir/ids_* | wc -l)

# If the user has an NCBI API key, we can send requests more frequently.
# See https://www.ncbi.nlm.nih.gov/books/NBK25497/ and man parallel.
#
# Temporarily turn off undefined variable usage errors because the user may
# not have the NCBI_API_KEY variable set.
set +u
if [ -z "$NCBI_API_KEY" ]
then
    delayArg='--delay 0.35s'
    apiParam=
else
    delayArg='--delay 0.1s'
    apiParam="api_key=${NCBI_API_KEY}&"
fi
set -u

failureCount=0
attempt=0
attemptLimit=5

# Note 5 repetitions of [a-z] here (see above comment).
lastFile=$(ls $idDir/ids_[a-z][a-z][a-z][a-z][a-z] | tail -n 1)

commandFile=$(mktemp)

if [ $dryRun -eq 0 ]
then
    # Clean up on exit.
    # Note 5 repetitions of [a-z] here (see above comment).
    trap 'rm -f $commandFile $idDir/ids_[a-z][a-z][a-z][a-z][a-z] $idDir/ids_[a-z][a-z][a-z][a-z][a-z].params' 0 1 2 3 15
fi

while [ $attempt -eq 0 -o $failureCount -gt 0 ]
do
    attempt=$((attempt + 1))
    
    if [ $attempt -gt $attemptLimit ]
    then
        echo "Attempt limit ($attemptLimit) reached. $failureCount GenBank files still not downloaded. Exiting." 2>&1
        break
    elif [ $attempt -gt 1 ]
    then
        echo "Re-running (attempt $attempt of $attemptLimit) to retry $failureCount incompletely downloaded files." 2>&1
    fi

    failureCount=0

    # Note 5 repetitions of [a-z] here (see above comment).
    for i in $idDir/ids_[a-z][a-z][a-z][a-z][a-z]
    do
        out=$i.gb

        if [ -f $out ]
        then
            # We know how many GenBank records there should be in each file.
            if [ $i = $lastFile ]
            then
                expectedCountArg=$(wc -l < $i | awk '{print $1}')
            else
                expectedCountArg=$maxIdsPerRequest
            fi

            if parse-genbank-flat-file.py --quiet --expectedCount $expectedCountArg $out
            then
                # $out already exists, can be parsed, and has the right number of records. Skip.
                continue
            else
                echo "$out already exists but cannot be parsed. Removing and re-downloading." >&2
                rm -f $out
                failureCount=$((failureCount + 1))
            fi
        fi

        echo -n "${apiParam}db=nucleotide&rettype=gb&id=" > $i.params
        cat $i | tr '\n' , | sed -e 's/,$//' >> $i.params
        echo "curl --silent -d @$i.params $URL > $out"
    done > $commandFile

    if [ $dryRun -eq 0 ]
    then
        parallel --bar $delayArg < $commandFile
    fi
done

if [ $dryRun -eq 1 ]
then
    echo "$count Genbank records would be downloaded"
    echo "The dry run id directory was $idDir"
    echo "Download commands are in $commandFile"
    # Just show what we'd do. We could instead cat the command file,
    # but that can be very verbose. Instead just let the user see the
    # file name and they can look in it if they care.
    echo "These would be executed using: parallel --bar $delayArg < $commandFile"
    echo "To clean up, run: rm -r '$idDir' '$commandFile'"
else
    echo "$count Genbank records downloaded"
fi

exit $failureCount

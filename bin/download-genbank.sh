#!/bin/bash

# This script downloads flat files (in GenBank format) of sequence info
# from NCBI.
#
# It is one of a set of utilities that we use to make a protein database
# based on the *clustered* nucleotide genomes found in the RVDB database
# found at https://hive.biochemistry.gwu.edu/rvdb
#
# Usage is e.g.:
#
#     download-genbank.sh [-n] ID_DIR C-RVDBv15.1.fasta ACCESSION-FILE
#
# where:
#
#     ID_DIR is the name of a directory in which to store files of
#         sequence ids, parameters to cURL, and GenBank files downloaded
#         from NCBI.
#
#     C-RVDBv15.1.fasta is the name of an RVDB clustered nucleotide FASTA file.
#
#     ACCESSION-FILE is the name of a file to write accession numbers and
#         sequence names to (this will also be needed by a later stage of
#         processing).
#
# The -n option can be used to make the script just print out what it would
# download, rather than actually doing the downloads. The ids directory and
# accession number to name file are still created.
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

# These are accession numbers in the RVDB database that correspond to dead
# (suppressed or withdrawn) sequence records (see ../doc/protein-database.md)
# for more info.
DEAD_REGEX='JX873962|KX912841|NC_021196|NC_036584'

while [ $# -gt 0 ]
do
    case "$1" in
        -n)
            dryRunArg=--dry-run
            shift
            ;;
        *)
            break
            ;;
    esac
done

case $# in
    3) idDir=$1; rvdb=$2; accessions=$3;
    *) echo "Usage: $(basename $0) id-dir RVDB-file.fasta accession-file" >&2; exit 1;;
esac

if [ ! -d $idDir ]
then
    mkdir "$idDir"
fi

trap 'rm -f $idDir/ids_[a-z][a-z][a-z][a-z][a-z] $idDir/ids_[a-z][a-z][a-z][a-z][a-z].params' 0 1 2 3 15

set +o pipefail
egrep '^>' $rvdb | cut -f3,4 -d\| | egrep -v "$DEAD_REGEX_REGEX" > $accessions
set -o pipefail

cut -f1 -d\| < $accessions | split -a 5 -l $maxIdsPerRequest - $idDir/ids_

# If the user has an NCBI API key, we can send requests more frequently.
# See https://www.ncbi.nlm.nih.gov/books/NBK25497/ and man parallel.
if [ -z "$NCBI_API_KEY" ]
then
    delayArg='--delay 0.35s'
    apiParam=
else
    delayArg='--delay 0.1s'
    apiParam="api_key=${NCBI_API_KEY}&"
fi

failureCount=0
attempt=0
attemptLimit=5

lastFile=$(ls $idDir/ids_[a-z][a-z][a-z][a-z][a-z] | tail -n 1)

while [ $attempt -eq 0 -o $failureCount -gt 0 ]
do
    attempt=$((attempt + 1))
    
    if [ $attempt -gt $attemptLimit ]
    then
        echo "Attempt limit ($attemptLimit) reached. $failureCount .genbank files still not downloaded. Exiting." 2>&1
        break
    elif [ $attempt -gt 1 ]
    then
        echo "Re-running (attempt $attempt of $attemptLimit) to retry $failureCount incompletely downloaded files." 2>&1
    fi

    failureCount=0

    for i in $idDir/ids_[a-z][a-z][a-z][a-z][a-z]
    do
        out=$i.genbank

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
    done | parallel --bar $dryRunArg $delayArg
done

exit $failureCount

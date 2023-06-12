#!/bin/sh

set -Eeuo pipefail

URL=https://ftp.ncbi.nlm.nih.gov/refseq/release/viral

foundSomething=0

# The following gets the list of files from $URL and pulls out the names we
# want. I do it this way because the number of refseq viral files does grow
# over time and I don't want to hard-code a maximum file number.
for file in $(curl -L --silent $URL | egrep '>viral\.[0-9]+\.genomic\.gbff\.gz' | cut -f2 -d\")
do
    echo "Downloading '$file'." >&2
    curl -L -O $URL/$file
    gunzip -t $file
    foundSomething=1
done

if [ $foundSomething -eq 0 ]
then
    echo "$(basename $0): Error! No refseq viral genome file names were identified at $URL." >&2
    echo "$(basename $0): Perhaps your version of egrep is failing to match the names." >&2
    exit 1
fi

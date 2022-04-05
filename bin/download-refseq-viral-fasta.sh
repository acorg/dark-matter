#!/bin/sh

set -Eeuo pipefail

URL=https://ftp.ncbi.nlm.nih.gov/refseq/release/viral

# The following gets the list of files from $URL and pulls out the names we
# want. I do it this way because the number of refseq viral files does grow
# over time and I don't want to hard-code a maximum file number.
for file in $(curl -L --silent $URL | egrep '>viral\.\d+\.\d+\.genomic\.fna\.gz' | cut -f2 -d\")
do
    echo "Downloading '$file'." >&2
    curl -L -O $URL/$file
    gunzip -t $file
done

cat viral.*.genomic.fna.gz | gunzip > viral.genomic.fasta
egrep '^>' < viral.genomic.fasta > viral.genomic.ids

#!/usr/bin/env python

import argparse
import requests

# See https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
# for URL format details.

URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    "db=%s&id=%s&rettype=%s&retmode=text"
)

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Fetch a sequence by id from NCBI and write it to stdout " "as FASTA."
    ),
)

parser.add_argument("id", help="The id of the sequence to fetch.")

parser.add_argument(
    "--database",
    default="nucleotide",
    choices=("nucleotide", "protein"),
    help="The name of the NCBI database to query.",
)

parser.add_argument(
    "--format", default="fasta", choices=("fasta", "gb"), help="The format to return."
)

args = parser.parse_args()

print(requests.get(URL % (args.database, args.id, args.format)).text.rstrip("\n"))

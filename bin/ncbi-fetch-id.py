#!/usr/bin/env python

from __future__ import print_function, division

import argparse
import requests

# See https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
# for URL format details.

URL = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
       'db=%(database)s&id=%(id)s&rettype=fasta&retmode=text')

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Fetch a sequence by id from NCBI and write it to stdout '
                 'as FASTA.'))

parser.add_argument(
    'id', help='The id of the sequence to fetch.')

parser.add_argument(
    '--database', default='nucleotide', choices=('nucleotide', 'protein'),
    help='The name of the NCBI database to query.')

args = parser.parse_args()

print(requests.get(
    URL % {'database': args.database, 'id': args.id}
).text.rstrip('\n'))

#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
import re

from dark.taxonomy import AccessionLineageFetcher, formatLineage

VERSION_REGEX = re.compile(r'\.\d+$')


def taxonomyInfo(accession, db, namesOnly, separator):
    """
    Convert the lineage information for C{accession} into a formatted string.

    @param accession: The C{str} accession number to look up.
    @param db: The sqlite3 taxonomy database to consult.
    @param namesOnly: If C{True} only print taxonomic names.
    @param separator: The C{str} separator to put between fields. If C{None}
        return space-padded aligned columns.
    @return: A formatted C{str} for printing or C{None} if C{accession} is
        not found in the taxonomy database.
    """
    lineage = db.lineage(accession)

    if lineage is None and VERSION_REGEX.search(accession) is None:
        # Try adding a version number.
        lineage = db.lineage(accession + '.1')

    if lineage:
        return formatLineage(lineage, namesOnly, separator)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Print taxonomy information for accession numbers')

    parser.add_argument(
        'accessions', nargs='*',
        help=('The accession number(s) to print taxonomy information for. If '
              'not given, accession numbers are read from standard input, one '
              'per line.'))

    parser.add_argument(
        '--separator',
        help=('The string to separate output columns with. Default is to '
              'print aligned output padded with multiple spaces (or a TAB if '
              '--namesOnly is used).'))

    parser.add_argument(
        '--database', required=True,
        help=('The file holding the sqlite3 taxonomy database. See '
              'https://github.com/acorg/ncbi-taxonomy-database for how to '
              'build one.'))

    parser.add_argument(
        '--namesOnly', default=False, action='store_true',
        help='If specified, only print the taxonomic names.')

    parser.add_argument(
        '--printAccession', default=False, action='store_true',
        help='If specified, also print the accession number.')

    args = parser.parse_args()

    db = AccessionLineageFetcher(args.database)

    if args.accessions:
        accessions = args.accessions
    else:
        accessions = (line[:-1] for line in sys.stdin)

    separator = args.separator

    if args.namesOnly and separator is None:
        separator = '\t'

    printAccession = args.printAccession

    for accession in accessions:
        if printAccession:
            print(accession + ':')
        result = taxonomyInfo(accession, db, args.namesOnly, separator)
        if result:
            print(result)
        else:
            print('%r was not found in the taxonomy database' % accession,
                  file=sys.stderr)

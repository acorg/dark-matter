#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
import re

from dark.taxonomy import (
    addTaxonomyDatabaseCommandLineOptions,
    parseTaxonomyDatabaseCommandLineOptions)

VERSION_REGEX = re.compile(r'\.\d+$')


def hosts(id_, db):
    """
    Look up hosts for an accession or taxonomy id.

    @param id_: The C{str} accession number or taxonomy id to look up.
    @param db: The sqlite3 taxonomy database to consult.
    @return: A C{set} of C{str} host names.
    """
    try:
        idInt = int(id_)
    except ValueError:
        hosts = db.hosts(id_)
    else:
        hosts = db.hosts(idInt)

    if hosts is None and VERSION_REGEX.search(id_) is None:
        # Try adding a version number.
        hosts = db.hosts(id_ + '.1')

    return hosts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Print hosts for accession numbers or taxonomy ids')

    parser.add_argument(
        'ids', nargs='*',
        help=('The ids (accession numbers, names, or taxonomy ids) to print '
              'host information for. If not given, ids are read from '
              'standard input, one per line.'))

    parser.add_argument(
        '--database', required=True,
        help=('The file holding the sqlite3 taxonomy database. See '
              'https://github.com/acorg/ncbi-taxonomy-database for how to '
              'build one.'))

    parser.add_argument(
        '--printId', default=False, action='store_true',
        help='If specified, also print the id.')

    addTaxonomyDatabaseCommandLineOptions(parser)

    args = parser.parse_args()

    db = parseTaxonomyDatabaseCommandLineOptions(args, parser)

    if args.ids:
        ids = args.ids
    else:
        ids = (line[:-1] for line in sys.stdin)

    for id_ in ids:
        if args.printId:
            print(id_ + ':')
        hosts = hosts(id_, db)
        if hosts:
            print(', '.join(sorted(hosts)))
        else:
            print('No host information for %r found in the taxonomy database.'
                  % id_, file=sys.stderr)

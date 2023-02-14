#!/usr/bin/env python

import sys
import argparse
import re

from dark.taxonomy import (
    formatLineage,
    addTaxonomyDatabaseCommandLineOptions,
    parseTaxonomyDatabaseCommandLineOptions,
)

VERSION_REGEX = re.compile(r"\.\d+$")


def taxonomyInfo(id_, db, namesOnly, separator):
    """
    Convert the lineage information for C{id_} into a formatted string.

    @param id_: The C{str} accession number or C{int} taxonomy id to
        look up.
    @param db: The sqlite3 taxonomy database to consult.
    @param namesOnly: If C{True} only print taxonomic names.
    @param separator: The C{str} separator to put between fields. If C{None}
        return space-padded aligned columns.
    @return: A formatted C{str} for printing or C{None} if C{id_} is
        not found in the taxonomy database.
    """
    try:
        idInt = int(id_)
    except ValueError:
        lineage = db.lineage(id_)
    else:
        lineage = db.lineage(idInt)

    if lineage is None and VERSION_REGEX.search(id_) is None:
        # Try adding a version number.
        lineage = db.lineage(id_ + ".1")

    if lineage:
        return formatLineage(lineage, namesOnly, separator)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Print taxonomy information for accession numbers or " "taxonomy ids."
        ),
    )

    parser.add_argument(
        "ids",
        nargs="*",
        metavar="id",
        help=(
            "The ids (accession numbers, names, or taxonomy ids) to print "
            "taxonomy information for. If not given, ids are read from "
            "standard input, one per line."
        ),
    )

    parser.add_argument(
        "--separator",
        help=(
            "The string to separate output columns with. Default is to "
            "print aligned output padded with multiple spaces (or a TAB if "
            "--namesOnly is used)."
        ),
    )

    parser.add_argument(
        "--namesOnly",
        default=False,
        action="store_true",
        help="If specified, only print the taxonomic names.",
    )

    parser.add_argument(
        "--print",
        default=False,
        action="store_true",
        help="If specified, also print the id.",
    )

    addTaxonomyDatabaseCommandLineOptions(parser)

    args = parser.parse_args()

    db = parseTaxonomyDatabaseCommandLineOptions(args, parser)

    if args.ids:
        ids = args.ids
    else:
        ids = (line.strip() for line in sys.stdin)

    separator = args.separator

    if args.namesOnly and separator is None:
        separator = "\t"

    for id_ in ids:
        if args.print:
            print(id_ + ":")
        result = taxonomyInfo(id_, db, args.namesOnly, separator)
        if result:
            print(result)
        else:
            print("%r was not found in the taxonomy database" % id_, file=sys.stderr)

#!/usr/bin/env python

import argparse

from dark.proteins import SqliteIndex


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Query an sqlite3 protein/genome database and write "
        "the matching sequences to standard output as FASTA."
    ),
)

parser.add_argument(
    "--databaseFile",
    required=True,
    help=(
        "The protein/genome database file. This must be an sqlite3 database "
        "file created by make-protein-database.py"
    ),
)

parser.add_argument(
    "--query",
    help=(
        "An SQL database query. "
        "The query must result in rows with two columns, the sequence id "
        "and the sequence. These can be protein or nucleotide sequences, "
        "but it is probably better not to mix them!"
    ),
)

parser.add_argument(
    "--queryFile",
    help=(
        "The name of a file containing an SQL database query. "
        "The query must result in rows with two columns, the sequence id "
        "and the sequence. These can be protein or nucleotide sequences, "
        "but it is probably better not to mix them!"
    ),
)

parser.add_argument(
    "--keepDuplicateIds",
    default=False,
    action="store_true",
    help="If specified, do not remove duplicate ids.",
)

parser.add_argument(
    "--idsOnly",
    default=False,
    action="store_true",
    help="If specified, only print matching sequence ids.",
)

args = parser.parse_args()

db = SqliteIndex(args.databaseFile)

if args.queryFile:
    query = open(args.queryFile).read()
elif args.query:
    query = args.query
else:
    raise RuntimeError("You must use either --query or --queryFile")

cur = db.execute(query)

removeDuplicates = not args.keepDuplicateIds
idsOnly = args.idsOnly

if removeDuplicates:
    seen = set()

while True:
    rows = cur.fetchmany()
    if rows:
        for row in rows:
            try:
                id_, sequence, genomeName = row
            except ValueError:
                id_, sequence = row
            else:
                id_ = "%s [%s]" % (id_, genomeName)

            if removeDuplicates:
                if id_ in seen:
                    continue
                else:
                    seen.add(id_)
            if idsOnly:
                print(id_)
            else:
                print(">%s\n%s" % (id_, sequence))
    else:
        break

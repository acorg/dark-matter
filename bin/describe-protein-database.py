#!/usr/bin/env python

import argparse

from dark.civ.proteins import SqliteIndex


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Print a description of a protein/genome database.",
)

parser.add_argument("databaseFile", help="The Sqlite3 database file.")

args = parser.parse_args()

db = SqliteIndex(args.databaseFile)

print("%d genomes with %d proteins" % (db.genomeCount(), db.proteinCount()))

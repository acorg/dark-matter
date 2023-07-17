#!/usr/bin/env python

import argparse
import bz2file
import sys

from dark.blast.conversion import XMLRecordsReader


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a BLAST XML file to JSON.",
        epilog=(
            "Give a BLAST XML file and convert it to JSON. Optionally "
            "compress the JSON output."
        ),
    )

    parser.add_argument(
        "--json",
        metavar="JSON-output-file",
        help=(
            "the JSON filename to write the converted XML to. If omitted, "
            "standard output will be used."
        ),
    )

    parser.add_argument(
        "--xml",
        metavar="BLAST-XML-file",
        default=sys.stdin,
        help=(
            "the BLAST XML output file to convert. If omitted, standard "
            "input will be read."
        ),
    )

    parser.add_argument(
        "--bzip2",
        default=False,
        action="store_true",
        help="If True, compress the JSON output using bzip2.",
    )

    args = parser.parse_args()

    if args.bzip2:
        fp = bz2file.BZ2File(args.json or sys.stdout, "w")
    else:
        fp = open(args.json, "w") if args.json else sys.stdout

    reader = XMLRecordsReader(args.xml)
    reader.saveAsJSON(fp)
    fp.close()

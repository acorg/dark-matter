#!/usr/bin/env python

import argparse

from dark.sam import samReferencesToStr

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Print SAM/BAM file reference names and lengths.",
)

parser.add_argument(
    "samfile", metavar="FILENAME", help="The name of a SAM/BAM alignment file."
)

print(samReferencesToStr(parser.parse_args().samfile))

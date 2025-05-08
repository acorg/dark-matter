#!/usr/bin/env python

"""
Read a file of our JSON BLAST output (which has one JSON object per
line) from stdin and pretty print it to stdout.
"""

import sys
from json import dumps, loads

for line in sys.stdin:
    s = dumps(loads(line[:-1]), sort_keys=True, indent=2)
    print("\n".join([field.rstrip() for field in s.splitlines()]))

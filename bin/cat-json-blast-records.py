#!/usr/bin/env python

"""
Read a file of our JSON BLAST output (which has one JSON object per
line) from stdin and pretty print it to stdout.
"""

from __future__ import print_function

from json import dumps, loads
import sys

for line in sys.stdin:
    s = dumps(loads(line[:-1]), sort_keys=True, indent=2)
    print('\n'.join([l.rstrip() for l in s.splitlines()]))

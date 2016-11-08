#!/usr/bin/env python

"""
Read summary output produced by group-summary-proteins.py from many samples
and group it by virus.

Reads filenames from standard input, writes to standard output.

The filenames on stdin will typically be the 05-panel/summary-virus files.
E.g.,

    $ find . -name summary-virus | group-summary-viruses.py
"""

from __future__ import print_function

import sys
import re
from collections import defaultdict

SUMMARY_RE = re.compile(' (\(\d+ proteins?, \d+ reads?\))$')

viruses = defaultdict(list)
fileCount = 0

for filename in sys.stdin:
    filename = filename[:-1]
    fileCount += 1
    with open(filename) as fp:
        for line in fp:
            line = line[:-1]
            match = SUMMARY_RE.search(line)
            if match:
                virus = line[:-len(match.group(0))]
                counts = match.group(1)

                simpleFilename = filename
                if simpleFilename.startswith('./'):
                    simpleFilename = simpleFilename[2:]
                if simpleFilename.endswith('/05-panel/summary-virus'):
                    simpleFilename = simpleFilename[:-23]

                viruses[virus].append('  %s %s' % (simpleFilename, counts))
            elif line:
                viruses[virus].append('    %s' % (line,))

print('Proteins from %d viruses were found across %d samples.' %
      (len(viruses), fileCount))
print()

for virus in sorted(viruses):
    fileInfo = viruses[virus]
    sampleCount = len([f for f in fileInfo if not f.startswith('    ')])
    print('%s (in %d sample%s)' % (virus, sampleCount,
                                   '' if len(fileInfo) == 1 else 's'))
    print('\n'.join(fileInfo))
    print()

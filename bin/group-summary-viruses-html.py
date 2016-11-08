#!/usr/bin/env python

"""
Read summary output produced by group-summary-proteins.py from many samples
and group it by virus.

Reads filenames from standard input, writes to standard output.

The filenames on stdin will typically be the 05-panel/summary-virus files.
E.g.,

    $ find . -name summary-virus | group-summary-viruses-html.py
"""

from __future__ import print_function

import sys
import re
from os.path import dirname, join
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

                viruses[virus].append(
                    '&nbsp;&nbsp;<tt>%s %s</tt> '
                    '(<a href="%s">panel</a>)' %
                    (simpleFilename, counts,
                     join(dirname(filename), 'out', 'index.html')))
            elif line:
                viruses[virus].append((line, filename))

print('<p>Proteins from %d viruses were found across %d samples.</p>' %
      (len(viruses), fileCount))
print()

for virus in sorted(viruses):
    samples = viruses[virus]
    sampleCount = len([s for s in samples if not isinstance(s, tuple)])
    print('<p><strong>%s (in %d sample%s)</strong><br>' %
          (virus, sampleCount, '' if len(samples) == 1 else 's'))
    for sample in samples:
        if isinstance(sample, tuple):
            line, filename = sample
            fields = line.split('\t')
            index = int(fields[5])
            print(
                '&nbsp;&nbsp;&nbsp;&nbsp;<tt>%s</tt> '
                '(<a href="%s">blue plot</a>, <a href="%s">fasta</a>)<br>' %
                ('\t'.join(fields[:5] + fields[6:]),
                 join(dirname(filename), 'out', '%d.png' % index),
                 join(dirname(filename), 'out', '%d.fasta' % index)))
        else:
            print('%s<br>' % sample)
    print('</p>')

#!/usr/bin/env python

"""
Read summary output produced by noninteractive-alignment-panel.py and group
matched subjects.

This is currently only useful when you are matching against a subject protein
database whose titles look e.g., like this:

gi|820945251|ref|YP_009137096.1| envelope glycoprotein H [Human herpesvirus 1]
gi|820945301|ref|YP_009137146.1| virion protein US10 [Human herpesvirus 1]
gi|820945229|ref|YP_009137074.1| ubiquitin E3 ligase ICP0 [Human herpesvirus 1]

In this case, those three matched subjects are proteins from the same virus.
This script will gather those matches under their common "Human herpesvirus 1"
title and print them together.

Reads from standard input, writes to standard output.
"""

from __future__ import print_function

import sys
import re
from collections import defaultdict
from operator import itemgetter

SUBJECT_RE = re.compile('^([^\[]+)\[([^\]]+)\]$')

subjectTitles = defaultdict(list)
proteinsWithNoSubject = []

for index, proteinLine in enumerate(sys.stdin):
    proteinLine = proteinLine[:-1]
    (coverage, medianScore, bestScore, readCount, hspCount, subjectLength,
     subjectTitle) = proteinLine.split('\t')

    match = SUBJECT_RE.match(subjectTitle)
    if match:
        proteinTitle = match.group(1).strip()
        subjectTitle = match.group(2)
        subjectTitles[subjectTitle].append({
            'coverage': float(coverage),
            'index': index,
            'medianScore': float(medianScore),
            'bestScore': float(bestScore),
            'readCount': int(readCount),
            'hspCount': int(hspCount),
            'proteinTitle': proteinTitle,
        })
    else:
        proteinsWithNoSubject.append(proteinLine)

titleGetter = itemgetter('proteinTitle')
readCountGetter = itemgetter('readCount')

for subjectTitle in sorted(subjectTitles):
    proteinMatches = subjectTitles[subjectTitle]
    proteinCount = len(proteinMatches)
    totalReads = sum(readCountGetter(p) for p in proteinMatches)
    print('%s (%d protein%s, %d read%s)' %
          (subjectTitle,
           proteinCount, '' if proteinCount == 1 else 's',
           totalReads, '' if totalReads == 1 else 's'))
    proteinMatches.sort(key=titleGetter)
    for proteinMatch in proteinMatches:
        print('  %(coverage)f\t%(medianScore)f\t%(bestScore)f\t'
              '%(readCount)d\t%(hspCount)d\t%(index)d\t%(proteinTitle)s' %
              proteinMatch)
    print()

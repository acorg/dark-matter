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

0.77 47.00 48.10 5  5  74  gi|101105594| ubiquitin [Brazilian marseillevirus]
0.31 42.70 48.10 47 47 630 gi|313768007| protein BpV1_008c [Bathycoccus virus]
0.21 42.00 48.10 42 42 687 gi|313768010| protein BpV1_011c [Bathycoccus virus]
0.77 46.60 48.10 5  5  74  gi|327409793| ubiquitin [Lausannevirus]
0.33 42.70 48.10 48 48 624 gi|472342805| protein 70 [Micromonas pusilla virus]

Input line fields are

        coverage
        median score
        best score
        read count
        HSP count
        index
        protein title

The 'index' is the 0-based index of all the proteins matched by the sample.
It is used to make a link to the blue plot and FASTA file (also produced by
noninteractive-alignment-panel.py) for the match.
"""

from __future__ import print_function
import argparse
import sys

from dark.proteins import ProteinGrouper


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Group proteins by the virus they're from")

    parser.add_argument(
        'filenames', nargs='*', help='Filenames to read input from')

    parser.add_argument(
        '--sampleNameRegex', default=None,
        help=('A regex to match the sample name in input filenames.'))

    parser.add_argument(
        '--html', default=False, action='store_true',
        help='If given, output HTML.')

    args = parser.parse_args()

    grouper = ProteinGrouper(sampleNameRegex=args.sampleNameRegex)

    if args.filenames:
        filenames = args.filenames
    else:
        filenames = (line[:-1] for line in sys.stdin)

    for filename in filenames:
        with open(filename) as fp:
            grouper.addFile(filename, fp)

    print(grouper.toHTML() if args.html else grouper.toStr())

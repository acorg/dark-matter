#!/usr/bin/env python

"""
Read protein match output produced by noninteractive-alignment-panel.py and
group it by virus.

This is currently only useful when you are matching against a subject protein
database whose titles have a virus name in square brackets, like this:

gi|820945251|ref|YP_009137096.1| envelope glycoprotein H [Human herpesvirus 1]
gi|820945301|ref|YP_009137146.1| virion protein US10 [Human herpesvirus 1]
gi|820945229|ref|YP_009137074.1| ubiquitin E3 ligase ICP0 [Human herpesvirus 1]

In this case, those three matched subjects are from the same virus. This script
will gather those matches under their common "Human herpesvirus 1" title and
provides methods to print them.

The script reads *file names* from standard input, and writes to standard
output.  Alternately, you can also provide file names on the command line.

Typical usage:

  $ find . -name summary-proteins | group-summary-proteins.py \
        --sampleNameRegex '(Sample_\d+)/' --html > index.html

Input files must contain lines in the following format:

0.77 47.00 48.10  5  5  74 gi|101105594| ubiquitin [Brazilian marseillevirus]
0.31 42.70 48.10 47 47 630 gi|313768007| protein BpV1_008c [Bathycoccus virus]
0.21 42.00 48.10 42 42 687 gi|313768010| protein BpV1_011c [Bathycoccus virus]
0.77 46.60 48.10  5  5  74 gi|327409793| ubiquitin [Lausannevirus]
0.33 42.70 48.10 48 48 624 gi|472342805| protein 70 [Micromonas pusilla virus]

Fields must be whitespace separated. The seven fields are:

    Coverage
    Median bit score
    Best bit score
    Read count
    HSP count
    Protein length
    Protein title (may contain whitespace)
"""

from __future__ import print_function
import argparse
import sys

# It's not clear that the PDF backend is the right choice here, but it
# works (i.e., the generation of PNG images works fine).
import matplotlib
matplotlib.use('PDF')

# These imports are here because dark.proteins imports matplotlib.pyplot
# and we need to set the matplotlib backend before the import. So please
# don't move this import higher in this file.

from dark.proteins import ProteinGrouper


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Group proteins by the virus they're from.")

    parser.add_argument(
        'filenames', nargs='*', help='Sample file names to read input from.')

    parser.add_argument(
        '--sampleNameRegex', default=None,
        help=('An (optional) regular expression that can be used to extract a '
              'short sample name from full sample file name.  The regular '
              'expression must have a matching group (delimited by '
              'parentheses) to capture the part of the file name that should '
              'be used as the sample name.'))

    parser.add_argument(
        '--virusPanelFilename', default=None,
        help=('An (optional) filename to write a virus-sample panel PNG '
              'image to.'))

    parser.add_argument(
        '--html', default=False, action='store_true',
        help='If specified, output HTML instead of plain text.')

    args = parser.parse_args()

    grouper = ProteinGrouper(sampleNameRegex=args.sampleNameRegex)

    if args.filenames:
        filenames = args.filenames
    else:
        filenames = (line[:-1] for line in sys.stdin)

    for filename in filenames:
        with open(filename) as fp:
            grouper.addFile(filename, fp)

    if args.html:
        print(grouper.toHTML(args.virusPanelFilename))
    else:
        print(grouper.toStr())

#!/usr/bin/env python

"""
Read protein match output produced by noninteractive-alignment-panel.py and
group it by pathogen (either virus or bacteria).

This is currently only useful when you are matching against a subject protein
database whose titles have a pathogen name in square brackets, like this:

gi|820945251|ref|YP_009137096.1| envelope glycoprotein H [Human herpesvirus 1]
gi|820945301|ref|YP_009137146.1| virion protein US10 [Human herpesvirus 1]
gi|820945229|ref|YP_009137074.1| ubiquitin E3 ligase ICP0 [Human herpesvirus 1]

In this case, those three matched subjects are from the same pathogen. This
script will gather those matches under their common "Human herpesvirus 1"
title and provides methods to print them.

Files with names in this format are provided by NCBI for their viral and
bacterial refseq protein FASTA.

The script reads *file names* from standard input, and writes to standard
output.  Alternately, you can also provide file names on the command line.

Typical usage:

  $ find . -name summary-proteins | proteins-to-pathogens.py \
        --sampleNameRegex "(Sample_\\d+)/" --html > index.html

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
    Title (in the format "protein name [pathogen name]")
"""

from __future__ import print_function
import argparse
import sys
from itertools import chain

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
        description="Group proteins by the pathogen they're from.")

    parser.add_argument(
        'filenames', nargs='*', help='Sample file names to read input from.')

    parser.add_argument(
        '--sampleName',
        help=('An (optional) sample name. This is only used in producing '
              'HTML output. Should be used when all input files are for a '
              'single sample.  Cannot be used with --sampleNameRegex.'))

    parser.add_argument(
        '--sampleNameRegex',
        help=('An (optional) regular expression that can be used to extract a '
              'short sample name from full sample file name.  The regular '
              'expression must have a matching group (delimited by '
              'parentheses) that captures the part of the file name that '
              'should be used as the sample name.'))

    parser.add_argument(
        '--pathogenPanelFilename',
        help=('An (optional) filename to write a pathogen-sample panel PNG '
              'image to.'))

    parser.add_argument(
        '--sampleIndexFilename',
        help=('An (optional) filename to write a sample index file to. '
              'Lines in the file will have an integer index, a space, and '
              'then the sample name. Only produced if --html is used '
              '(because the pathogen-NNN-sample-MMM.fastq are only written '
              'in that case).'))

    parser.add_argument(
        '--pathogenIndexFilename',
        help=('An (optional) filename to write a pathogen index file to. '
              'Lines in the file will have an integer index, a space, and '
              'then the pathogen name. Only produced if --html is used '
              '(because the pathogen-NNN-sample-MMM.fastq are only written '
              'in that case).'))

    parser.add_argument(
        '--html', default=False, action='store_true',
        help='If specified, output HTML instead of plain text.')

    parser.add_argument(
        '--format', default='fasta', choices=('fasta', 'fastq'),
        help=('Give the format of the sequence files written by '
              'noninteractive-alignment-panel.py when it created the '
              'summary-proteins files given on output.'))

    parser.add_argument(
        '--proteinFastaFilename', '--pff', nargs='+', action='append',
        help=('Optional filename(s) giving the name of the FASTA file(s) '
              'with the protein AA sequences with their associated pathogens '
              'in square brackets. This is the format used by NCBI for '
              'bacterial and viral reference sequence protein files. If '
              'given, the contents of this file will be used to determine how '
              'many proteins each matched pathogen has. This makes it much '
              'easier to spot significant matches (as opposed to those where, '
              'say, just one protein from a pathogen is matched).'))

    parser.add_argument(
        '--minProteinFraction', type=float, default=0.0,
        help=('The minimum fraction of proteins in a pathogen that must be '
              'matched by a particular sample in order for that pathogen to '
              'be displayed for that sample.'))

    parser.add_argument(
        '--pathogenType', default='viral', choices=('bacterial', 'viral'),
        help=('Specify the pathogen type. This option only affects the '
              'language used in HTML output.'))

    parser.add_argument(
        '--showReadLengths', default=False, action='store_true',
        help=('If specified, the HTML output (use --html to get this) will '
              'contain the lengths of all reads that match proteins for a '
              'pathogen.'))

    parser.add_argument(
        '--assetDir', default='out',
        help=('The output directory where noninteractive-alignment-panel.py '
              'puts its HTML, plots and FASTA or FASTQ files, needed for '
              'using --html'))

    args = parser.parse_args()

    if args.sampleName and args.sampleNameRegex:
        print('It does not make sense to use --sampleName '
              'as well as --sampleNameRegex', file=sys.stderr)
        sys.exit(1)

    if not args.html:
        if args.sampleIndexFilename:
            print('It does not make sense to use --sampleIndexFilename '
                  'without also using --html', file=sys.stderr)
            sys.exit(1)
        if args.pathogenIndexFilename:
            print('It does not make sense to use --pathogenIndexFilename '
                  'without also using --html', file=sys.stderr)
            sys.exit(1)

    if args.proteinFastaFilename:
        # Flatten lists of lists that we get from using both nargs='+' and
        # action='append'. We use both because it allows people to use
        # (e.g.)  --pff on the command line either via "--pff file1 --pff
        # file2" or "--pff file1 file2", or a combination of these. That
        # way it's not necessary to remember which way you're supposed to
        # use it and you also can't be hit by the subtle problem
        # encountered in https://github.com/acorg/dark-matter/issues/453
        proteinFastaFilenames = list(chain.from_iterable(
            args.proteinFastaFilename))
    else:
        proteinFastaFilenames = None

    grouper = ProteinGrouper(assetDir=args.assetDir,
                             sampleName=args.sampleName,
                             sampleNameRegex=args.sampleNameRegex,
                             format_=args.format,
                             proteinFastaFilenames=proteinFastaFilenames,
                             saveReadLengths=args.showReadLengths)

    if args.filenames:
        filenames = args.filenames
    else:
        filenames = (line[:-1] for line in sys.stdin)

    for filename in filenames:
        with open(filename) as fp:
            grouper.addFile(filename, fp)

    if args.html:
        print(grouper.toHTML(args.pathogenPanelFilename,
                             minProteinFraction=args.minProteinFraction,
                             pathogenType=args.pathogenType,
                             sampleIndexFilename=args.sampleIndexFilename,
                             pathogenIndexFilename=args.pathogenIndexFilename))
    else:
        print(grouper.toStr())

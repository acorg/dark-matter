#!/usr/bin/env python

# See https://samtools.github.io/hts-specs/SAMv1.pdf for the SAM file
# format specification.

from __future__ import print_function, division

import sys
import argparse
from os.path import join
from resource import getrlimit, RLIMIT_NOFILE
from functools import partial

from dark.diamond.conversion import FIELDS
from dark.diamond.sam import (
    SimpleDiamondSAMWriter, PerReferenceDiamondSAMWriter)
from dark.proteins import SqliteIndex

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Convert DIAMOND tabular format to SAM. The DIAMOND '
                 'invocation *must* include --outfmt 6 %s' % FIELDS))

parser.add_argument(
    '--printFields', default=False, action='store_true',
    help=('Print the field names in the order that they must be given to '
          'diamond --outfmt 6 to produce correct input for this script, '
          'then exit.'))

parser.add_argument(
    '--mappingQuality', type=int, default=255,
    help=('The mapping quality to use for MAPQ (field 5). The default (255) '
          'indicates that mapping quality information is not available.'))

parser.add_argument(
    '--ram', action='store_true', default=False,
    help=('Do not use a temporary file to hold the non-header SAM output. '
          'This will run faster but use more memory since all non-header SAM '
          'output will be stored in RAM and only written out when the full '
          'header can be determined.'))

parser.add_argument(
    '--keepDescriptions', action='store_true', default=False,
    help=('Do not discard text after the first space in query or subject '
          'sequence ids. Note that this violates the SAM specification, but '
          'since SAM files are TAB-separated it may not cause problems '
          'downstream, depending on what tools you use on the SAM.'))

parser.add_argument(
    '--genomesProteinsDatabaseFilename', required=True,
    help=('The name of the file containing an sqlite3 database with '
          'information on genomes and proteins, as made by '
          'make-protein-database.py'))

parser.add_argument(
    '--outputDirectory', default='.',
    help=('The name of the directory to create per-reference SAM files in. '
          'This is only used when --perReferenceOutput is given, otherwise '
          'the output is written to stdout as a single SAM file.'))

parser.add_argument(
    '--perReferenceOutput', action='store_true', default=False,
    help=('If given, a separate SAM file is written for matches against each '
          'reference sequence. The files will have names corresponding to '
          'the accession number of the reference, with a .sam suffix '
          'and will be put into the directory specified by '
          '--outputDirectory.'))

parser.add_argument(
    '--filePointerCacheMaxSize', type=int,
    default=getrlimit(RLIMIT_NOFILE)[0] >> 1,
    help=('The maximum size for the file pointer cache. SAM output may '
          'be being written to an arbitrary number of output files so we '
          'need to make sure the OS limit on open files is not exceeded. '
          'The default value is size half the maximum number of open files '
          'will be used.'))

args = parser.parse_args()

if args.printFields:
    print(FIELDS)
    sys.exit(0)

if 0 > args.mappingQuality > 255:
    raise ValueError('Mapping quality must be between 0 and 255 (inclusive)')

genomesProteins = SqliteIndex(args.genomesProteinsDatabaseFilename)

if args.perReferenceOutput:

    def referenceTitleToBasename(subjectTitle):
        """
        Make a basename for a SAM output file for a DIAMOND subject match.

        @param subjectTitle: A C{str} database subject id. This must have a
            form like "acc|GENBANK|YP_009137150.1|GENBANK|NC_001798.2|..."
            where the accession number of the matched sequence is in the third
            '|'-separated field (which must not contain any characters that
            cannot be used in a filename).
        @return: A C{str} pathname to store the SAM for the reference to.
        """
        return join(args.outputDirectory,
                    genomesProteins.genomeAccession(subjectTitle))

    writer = PerReferenceDiamondSAMWriter(
        genomesProteins, mappingQuality=args.mappingQuality,
        ram=args.ram, keepDescriptions=args.keepDescriptions,
        baseFilenameFunc=referenceTitleToBasename,
        fpcMaxsize=args.filePointerCacheMaxSize)

    save = writer.save
else:
    writer = SimpleDiamondSAMWriter(
        genomesProteins, mappingQuality=args.mappingQuality,
        ram=args.ram, keepDescriptions=args.keepDescriptions)

    save = partial(writer.save, filename=sys.stdout)

addMatch = writer.addMatch

for line in sys.stdin:
    addMatch(line)

save()

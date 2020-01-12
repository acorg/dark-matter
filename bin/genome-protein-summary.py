#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from itertools import chain

from dark.civ.proteins import SqliteIndex
from dark.errors import NoSuchGenomeError
from dark.filter import (
    addFASTAFilteringCommandLineOptions, parseFASTAFilteringCommandLineOptions)
from dark.genbank import GenomeRanges
from dark.genomes import GenomeProteinInfo
from dark.reads import Reads
from dark.sam import SAMFilter, samReferences


def pct(a, b):
    assert a <= b
    if b:
        return ('%d/%d (%.3f%%)' %
                (a, b, (a / b if b else 0.0) * 100.0))
    else:
        return '0/0 (0.00%)'


def summarize(gpi, sortOn):
    """
    Print a summary of the genome proteins.

    @param gpi: A C{GenomeProteinInfo} instance.
    """
    genome = gpi.genome

    print('Summary of %s (%s):' % (genome['name'], genome['accession']))
    print('  Length: %d' % genome['length'])
    print('  Protein count: %d' % genome['proteinCount'])
    print('  Total protein offsets: %s' % (
        pct(len(gpi.offsets), genome['length'])))

    if gpi.samFiles:
        print('  SAM files analyzed: %d' % len(gpi.samFiles))
        for i, filename in enumerate(gpi.samFiles, start=1):
            print('    %d: %s' % (i, filename))
    else:
        return

    print('  Whole genome coverage (not just proteins):')
    print('    Reads matching genome: %d' % len(gpi.readIdsMatchingGenome))
    print('    Covered offsets: %s' % (
        pct(len(gpi.coveredOffsetCount), genome['length'])))
    print('    Average depth: %.3f' % (
        sum(gpi.coveredOffsetCount.values()) / genome['length']))

    coveredProteinOffsetCount = coveredProteinBasesCount = 0
    for offset in gpi.offsets:
        coveredProteinOffsetCount += bool(gpi.coveredOffsetCount[offset])
        coveredProteinBasesCount += gpi.coveredOffsetCount[offset]

    print('  Total protein coverage:')
    print('    Reads matching proteins: %d' % len(gpi.readIdsForAllProteins()))
    print('    Covered offsets: %s' % (
        pct(coveredProteinOffsetCount, len(gpi.offsets))))
    print('    Average depth: %.3f' % (
        coveredProteinBasesCount / len(gpi.offsets)))

    print('  Proteins matched: %s (sorted by %s):' % (
        pct(len(gpi.coveredProteins), genome['proteinCount']), sortOn))

    if sortOn == 'name':
        def key(proteinAccession):
            return gpi.proteins[proteinAccession]['name']
        reverse = False
    elif sortOn == 'offset':
        def key(proteinAccession):
            return GenomeRanges(
                gpi.proteins[proteinAccession]['offsets']).ranges[0][0]
        reverse = False
    elif sortOn == 'readCount':
        def key(proteinAccession):
            coverage = gpi.proteinCoverageInfo(proteinAccession)
            return len(coverage['readIds'])
        reverse = True
    elif sortOn == 'coverage':
        def key(proteinAccession):
            coverage = gpi.proteinCoverageInfo(proteinAccession)
            return coverage['coveredOffsets'] / coverage['ntLength']
        reverse = True
    elif sortOn == 'depth':
        def key(proteinAccession):
            coverage = gpi.proteinCoverageInfo(proteinAccession)
            return coverage['totalBases'] / coverage['ntLength']
        reverse = True

    for i, proteinAccession in enumerate(
            sorted(gpi.coveredProteins, key=key, reverse=reverse),
            start=1):
        protein = gpi.proteins[proteinAccession]
        print('    %d: %s (%d AA, %d nt with stop codon, %s)' %
              (i, protein['product'], protein['length'],
               protein['length'] * 3 + 3, protein['accession']))

        coverage = gpi.proteinCoverageInfo(proteinAccession)

        print('      Read count: %d' % len(coverage['readIds']))

        print('      Covered offsets: %s' % (
            pct(coverage['coveredOffsets'], coverage['ntLength'])))

        print('      Average depth: %.3f' % (
            coverage['totalBases'] / coverage['ntLength']))

        print('      Offsets: %s' % protein['offsets'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Print SAM/BAM file protein match statistics.')

    parser.add_argument(
        '--proteinGenomeDatabase', required=True,
        help=('The filename of an Sqlite3 database holding protein and '
              'genome information, as built by make-protein-database.py'))

    parser.add_argument(
        '--sortOn', default='readCount',
        choices=('coverage', 'depth', 'name', 'offset', 'readCount'),
        help='How to sort proteins for output.')

    parser.add_argument(
        '--skipTranslationChecks', dest='checkTranslations',
        action='store_false', default=True,
        help=('Skip the sanity check that database protein sequences can all '
              'be translated from the database genome sequence.'))

    addFASTAFilteringCommandLineOptions(parser)
    SAMFilter.addFilteringOptions(
        parser, samfileIsPositional=False, samfileAction='append',
        samfileNargs='*', samfileRequired=False)

    args = parser.parse_args()

    samfiles = list(chain.from_iterable(args.samfile)) if args.samfile else []

    if samfiles:
        if args.referenceId:
            referenceIds = args.referenceId
        else:
            # If all SAM files have just one reference and they're all the
            # same, use that. Else complain.
            referenceIds = set()
            for filename in samfiles:
                referenceIds.update(samReferences(filename))

            if len(referenceIds) != 1:
                print('No reference id(s) specified with --referenceId, and '
                      'the given SAM/BAM files do not contain exactly one '
                      '(identical) reference. Please use --referenceId')
                sys.exit(1)
    else:
        if args.referenceId:
            referenceIds = args.referenceId
        else:
            print('No reference id(s) specified with --referenceId.')
            sys.exit(1)

    # We don't have a file of reads, we just want a read filter that we
    # can use to filter the SAM file query sequences.
    reads = parseFASTAFilteringCommandLineOptions(args, Reads())
    samFilter = SAMFilter.parseFilteringOptions(args, reads.filterRead)
    filterAlignment = samFilter.filterAlignment

    proteinGenomeDB = SqliteIndex(args.proteinGenomeDatabase)

    for referenceId in referenceIds:
        try:
            gpInfo = GenomeProteinInfo(
                referenceId, proteinGenomeDB,
                checkTranslations=args.checkTranslations)
        except NoSuchGenomeError:
            print('Reference %r not found in genome database. Ignoring.' %
                  referenceId, file=sys.stderr)
        else:
            if samfiles:
                print('Processing %d SAM file%s for matches with %r:' %
                      (len(samfiles), '' if len(samfiles) == 1 else 's',
                       referenceId), file=sys.stderr)
                for i, filename in enumerate(samfiles, start=1):
                    print('  %d: %s' % (i, filename), file=sys.stderr)
                    gpInfo.addSAM(filename, filterAlignment)

            summarize(gpInfo, args.sortOn)

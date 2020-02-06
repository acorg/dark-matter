#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse
import requests
import multiprocessing
from os.path import basename, exists, join
from tempfile import mkdtemp

from dark.process import Executor


def saveStdin():
    if os.isatty(0):
        print('Reading sequences to match against from stdin.',
              file=sys.stderr)
    dirname = mkdtemp()
    filename = join(dirname, 'stdin.fastq')
    count = 0

    with open(filename, 'w') as fp:
        for line in sys.stdin:
            fp.write(line)
            count += 1

    if count == 0:
        print('Standard input was empty! Exiting.', file=sys.stderr)
        sys.exit(1)

    return dirname, filename


def makeIndexFromFastaFile(fastaFilename, tempdir, base, e):
    dirname = tempdir or mkdtemp()
    index = join(dirname, basename(base))

    e.execute("bowtie2-build --quiet '%s' '%s'" %
              (fastaFilename, index))

    return dirname, index


def makeIndexFromAccession(accessionId, tempdir, base, e, verbose, dryRun):
    dirname = tempdir or mkdtemp()
    fastaFilename = join(dirname, accessionId + '.fasta')
    index = join(dirname, basename(base))

    if not dryRun:
        if verbose:
            print('Downloading FASTA for accession %s from NCBI.' %
                  accessionId, file=sys.stderr)

        URL = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
               'db=%(database)s&id=%(id)s&rettype=fasta&retmode=text')

        with open(fastaFilename, 'w') as fp:
            print(requests.get(
                URL % {'database': 'nucleotide', 'id': accessionId}
            ).text.rstrip('\n'), file=fp)

    e.execute("bowtie2-build --quiet '%s' '%s'" %
              (fastaFilename, index))

    return dirname, index


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Run bowtie2 on a FASTA file. Optionally convert the result '
                 'to BAM, sorting, and indexing.'))

parser.add_argument(
    '--index', required=True,
    help=('Either: an accession number, a filename (ending in .fasta or .fa), '
          'or the name of a pre-existing bowtie2 index (created with '
          'bowtie2-build index).'))

parser.add_argument(
    '--query',
    help=('The FASTQ reads to match against the bowtie2 index given by '
          '--index. If not given, standard input will be read.'))

parser.add_argument(
    '--bowtie2Args', default='--end-to-end --no-unal',
    help=('Extra arguments to be passed to Bowtie2. Use --threads to specify '
          'a thread count.'))

parser.add_argument(
    # 3332 = UNMAP,SECONDARY,DUP,SUPPLEMENTARY
    '--samtoolsViewArgs', default='-F 3332 -q 30',
    help='Arguments to be passed to samtools view when creating the BAM file.')

parser.add_argument(
    '--base',
    help=('The base name of files to create. Suffixes such as .sam and .bam '
          'will be added. If not given, the prefix of the --query file '
          '(if any) will be used, stripped of its final suffix. If not given '
          'and no filename is given via --query, an error will be thrown.'))

parser.add_argument(
    '--picard',
    help=('The path to the Picard jar file. If given, Picard will be used to '
          'mark duplicates, which will then be removed using samtools. See '
          'https://github.com/broadinstitute/picard for details on Picard.'))

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help='Print the commands that were (or would be) executed.')

parser.add_argument(
    '--log', default=False, action='store_true',
    help='Show a log of commands that were (or would be) executed.')

parser.add_argument(
    '--threads', type=int,
    help='The number of threads to use when running bowtie2')

parser.add_argument(
    '--noBAM', default=False, action='store_true',
    help='If specified, do not convert SAM to BAM.')

parser.add_argument(
    '--noSort', default=False, action='store_true',
    help='If specified, do not sort the BAM.')

parser.add_argument(
    '--noIndex', default=False, action='store_true',
    help='If specified, do not index the sorted BAM.')

parser.add_argument(
    '--noClean', default=False, action='store_true',
    help=('If specified, do not remove intermediate .sam and .bam files or '
          'the temporary directory (if one is actually made).'))

parser.add_argument(
    '--force', default=False, action='store_true',
    help='If specified, overwrite pre-existing output files.')

parser.add_argument(
    '--dryRun', default=False, action='store_true',
    help='Do not run commands, just print what would be done.')

args = parser.parse_args()

if args.threads is None:
    nThreads = multiprocessing.cpu_count()
else:
    nThreads = args.threads


if args.query:
    tempDir, queryFile = None, args.query
else:
    tempDir, queryFile = saveStdin()


if args.base is None:
    if args.query is None:
        print('You must either use --base to specify an output file basename '
              'or else use --query to give an input name (whose suffix will '
              'be stripped and used as the output name.', file=sys.stderr)
        sys.exit(1)
    else:
        fields = args.query.rsplit('.', 1)
        if len(fields) < 2:
            print('No --base argument was given and the --query argument '
                  'does not have a .suffix that can be stripped.',
                  file=sys.stderr)
            sys.exit(1)
        base = fields[0]
else:
    base = args.base

e = Executor(args.dryRun)

# Find or make the bowtie2 index.

if exists(args.index):
    # Check if this is a pre-existing bowtie2 index. Look for and remove a
    # bowtie2 suffix, if any (otherwise bowtie2 will complain). We do
    # things this way to allow the user to use --index on the command line
    # and let TAB completion give them the full path to any bowtie index
    # file.
    for suffix in '1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2'.split():
        suffix = '.' + suffix
        if args.index.endswith(suffix):
            indexName = args.index[:-len(suffix)]
            break
    else:
        # Assume a FASTA file and make an index.
        tempDir, indexName = makeIndexFromFastaFile(
            args.index, tempDir, base, e)
else:
    # Not a filename. So either the start of the path to a bowtie2 index or
    # else an accession number.
    if exists(args.index + '.1.bt2'):
        indexName = args.index
    else:
        # Assume an accession number.
        tempDir, indexName = makeIndexFromAccession(
            args.index, tempDir, base, e, args.verbose, args.dryRun)

samFile = base + '.sam'
bamFile = base + '.bam'
sortedBamFile = base + '-sorted.bam'

if not (args.force or args.dryRun):
    existing = []
    for filename in samFile, bamFile, sortedBamFile:
        if exists(filename):
            existing.append(filename)
    if existing:
        print('Will not overwrite pre-existing file%s %s. '
              'Use --force to make me.' % (
                  '' if len(existing) == 1 else 's',
                  ', '.join(existing)),
              file=sys.stderr)
        sys.exit(2)

if args.verbose and not args.dryRun:
    print('Running bowtie2.', file=sys.stderr)

e.execute("bowtie2 %s --threads %d -x '%s' '%s' > '%s'" % (
    args.bowtie2Args, nThreads, indexName, queryFile, samFile))


if not args.noBAM:
    if args.verbose and not args.dryRun:
        print('Converting SAM to BAM.', file=sys.stderr)

    e.execute("samtools view -b %s < '%s' > '%s'" %
              (args.samtoolsViewArgs, samFile, bamFile))

    if not args.noClean:
        e.execute("rm '%s'" % samFile)

    if not args.noSort:
        if args.verbose and not args.dryRun:
            print('Sorting BAM.', file=sys.stderr)

        e.execute("samtools sort '%s' > '%s'" % (bamFile, sortedBamFile))

        e.execute("mv '%s' '%s'" % (sortedBamFile, bamFile))

        if args.picard:
            if args.verbose and not args.dryRun:
                print('Marking duplicates with Picard.', file=sys.stderr)

            tempDir = tempDir or mkdtemp()
            tempBAMFile = join(tempDir, 'picard.bam')
            tempErrFile = join(tempDir, 'picard.errs')
            e.execute('java -Xmn2g -Xms2g -Xmx2g -jar %s '
                      'MarkDuplicates I="%s" O="%s" M=/dev/null >"%s" 2>&1' %
                      (args.picard, bamFile, tempBAMFile, tempErrFile))

            if args.verbose and not args.dryRun:
                print('Removing duplicates marked by Picard.', file=sys.stderr)

            e.execute("samtools view -b -F 1024 < '%s' > '%s'" %
                      (tempBAMFile, bamFile))

        if not args.noIndex:
            if args.verbose and not args.dryRun:
                print('Indexing BAM.', file=sys.stderr)

            e.execute("samtools index '%s'" % bamFile)

if tempDir:
    if args.noClean:
        if args.verbose:
            print('Temporary directory was %r' % tempDir, file=sys.stderr)
    else:
        e.execute("rm -r '%s'" % tempDir)

if args.dryRun or args.log:
    print('\n'.join(e.log), file=sys.stderr)

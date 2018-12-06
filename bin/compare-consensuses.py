#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from os.path import exists, join
from os import mkdir
from tempfile import mkdtemp

from dark.process import Executor


def makeOuputDir(outputDir, force):
    """
    Create or check for an output directory.

    @param outputDir: A C{str} output directory name, or C{None}.
    @param force: If C{True}, allow overwriting of pre-existing files.
    @return: The C{str} output directory name.
    """
    if outputDir:
        if exists(outputDir):
            if not force:
                print('Will not overwrite pre-existing files. Use --force to '
                      'make me.', file=sys.stderr)
                sys.exit(1)
        else:
            mkdir(outputDir)
    else:
        outputDir = mkdtemp()
        print('Writing output files to %s' % outputDir)

    return outputDir


def samtoolsMpileup(outFile, referenceFile, alignmentFile, executor):
    """
    Use samtools mpileup to generate VCF.

    @param outFile: The C{str} name to write the output to.
    @param referenceFile: The C{str} name of the FASTA file with the reference
        sequence.
    @param alignmentFile: The C{str} name of the SAM or BAM alignment file.
    @param executor: An C{Executor} instance.
    """
    executor.execute(
        'samtools mpileup -u -v -f %s %s > %s' %
        (referenceFile, alignmentFile, outFile))


def bcftoolsMpileup(outFile, referenceFile, alignmentFile, executor):
    """
    Use bcftools mpileup to generate VCF.

    @param outFile: The C{str} name to write the output to.
    @param referenceFile: The C{str} name of the FASTA file with the reference
        sequence.
    @param alignmentFile: The C{str} name of the SAM or BAM alignment file.
    @param executor: An C{Executor} instance.
    """
    executor.execute(
        'bcftools mpileup -Ov -f %s %s > %s' %
        (referenceFile, alignmentFile, outFile))


def bcftoolsCallMulti(outFile, vcfFile, executor):
    """
    Use bcftools to make consensus calls, using the  -m (multiallelic) option.

    @param outFile: The C{str} name to write the output to.
    @param vcfFile: The C{str} name of the VCF file with the calls from
        the pileup.
    @param executor: An C{Executor} instance.
    """
    # Note that this just gives a call on the first genome in the calls file.
    executor.execute(
        'bcftools call -m -Ov -o %s < %s' % (outFile, vcfFile))


def bcftoolsCallConsensus(outFile, vcfFile, executor):
    """
    Use bcftools to make consensus calls, using the original -c option.

    @param outFile: The C{str} name to write the output to.
    @param vcfFile: The C{str} name of the VCF file with the calls from
        the pileup.
    @param executor: An C{Executor} instance.
    """
    # Note that this just gives a call on the first genome in the calls file.
    executor.execute(
        'bcftools call -c -Ov -o %s < %s' % (outFile, vcfFile))


def bcftoolsConsensus(outFile, vcfFile, id_, referenceFile, executor):
    """
    Use bcftools to extract consensus FASTA.

    @param outFile: The C{str} name to write the output to.
    @param vcfFile: The C{str} name of the VCF file with the calls from
        the pileup.
    @param id_: The C{str} identifier to use in the resulting FASTA sequence.
    @param referenceFile: The C{str} name of the FASTA file with the reference
        sequence.
    @param executor: An C{Executor} instance.
    """
    bgz = vcfFile + '.gz'
    executor.execute('bgzip -c %s > %s' % (vcfFile, bgz))
    executor.execute('tabix %s' % bgz)
    executor.execute(
        'bcftools consensus %s < %s | '
        'filter-fasta.py --idLambda \'lambda id: "%s"\' > %s' %
        (bgz, referenceFile, id_, outFile))


def vcfutilsConsensus(outFile, vcfFile, id_, _, executor):
    """
    Use vcftools to extract consensus FASTA.

    @param outFile: The C{str} name to write the output to.
    @param vcfFile: The C{str} name of the VCF file with the calls from
        the pileup.
    @param id_: The C{str} identifier to use in the resulting FASTA sequence.
    @param executor: An C{Executor} instance.
    """
    executor.execute(
        'vcfutils.pl vcf2fq < %s | '
        'filter-fasta.py --fastq --quiet --saveAs fasta '
        '--idLambda \'lambda id: "%s"\' > %s' %
        (vcfFile, id_, outFile))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compare different consensus making methods.')

    parser.add_argument(
        '--referenceFile', required=True, metavar='FILENAME',
        help='The name of the FASTA file containing the reference sequence.')

    parser.add_argument(
        '--alignmentFile', required=True, metavar='FILENAME',
        help=('The name of the SAM or BAM file containing an alignment to the '
              'reference.'))

    parser.add_argument(
        '--verbose', type=int, default=0, metavar='N',
        help=('The integer verbosity level (0 = no output, 1 = some output, '
              '2 = maximal output).'))

    parser.add_argument(
        '--force', default=False, action='store_true',
        help='If given, overwrite pre-existing files.')

    parser.add_argument(
        '--outputDir', metavar='DIRNAME',
        help='The directory to save result files to.')

    args = parser.parse_args()

    outputDir = makeOuputDir(args.outputDir, args.force)

    executor = Executor()

    pileuppers = (
        ('samtools-mpileup', samtoolsMpileup),
        ('bcftools-mpileup', bcftoolsMpileup))

    callers = (
        ('bcftools-c', bcftoolsCallConsensus),
        ('bcftools-m', bcftoolsCallMulti))

    consensusers = (
        ('bcftools-consensus', bcftoolsConsensus),
        ('vcfutils-vcf2fq', vcfutilsConsensus))

    consensusFiles = []

    for pileupName, pileupFunc in pileuppers:
        pileupMiddle = pileupName
        pileupFile = join(outputDir, 'pileup-' + pileupMiddle + '.vcf')
        pileupFunc(pileupFile, args.referenceFile,
                   args.alignmentFile, executor)

        for callerName, callerFunc in callers:
            callerMiddle = pileupMiddle + '-' + callerName
            callerFile = join(outputDir, 'calls-' + callerMiddle + '.vcf')
            callerFunc(callerFile, pileupFile, executor)

            for consensusName, consensusFunc in consensusers:
                consensusMiddle = callerMiddle + '-' + consensusName
                consensusFile = join(outputDir,
                                     'consensus-' + consensusMiddle + '.fasta')
                consensusFiles.append(consensusFile)
                consensusFunc(consensusFile, callerFile, consensusMiddle,
                              args.referenceFile, executor)

    # Let's assume there's at least one consensus file.
    consensusesFile = join(outputDir, 'consensuses.fasta')
    executor.execute(
        'cat %s > %s' % (' '.join(consensusFiles), consensusesFile))

    htmlFile = join(outputDir, 'consensus-identity.html')
    executor.execute(
        ('fasta-identity-table.py --footer --showGaps --showLengths < %s | '
         "perl -pe 's/-(bcftools|vcfutils)/ $1/g' > %s") %
        (consensusesFile, htmlFile))

    verbose = args.verbose
    if verbose > 0:
        print('The following commands were executed:')
        for line in executor.log:
            if line.startswith('#'):
                if verbose > 1:
                    print(line)
            else:
                print(line)

    print('Identity table comparing %d consensuses written to %s' %
          (len(consensusFiles), htmlFile))

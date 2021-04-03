#!/usr/bin/env python

import os
import sys
import argparse
from tempfile import mkdtemp
from os.path import join, basename

from dark.fasta import FastaReads
from dark.process import Executor

IVAR_FREQUENCY_THRESHOLD_DEFAULT = 0.6
IVAR_DOCS = (
    'https://andersen-lab.github.io/ivar/html/manualpage.html#autotoc_md19')


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Make a consensus sequence.')

    parser.add_argument(
        '--reference', required=True,
        help='The reference FASTA file.')

    parser.add_argument(
        '--bam',
        help=('The BAM file from which the consensus should be made. '
              'Required if --maskLowCoverage is used. If no BAM file is '
              'given, a VCF file must be provided. If both a BAM and a VCF '
              'file are given, the VCF file will take precedence.'))

    parser.add_argument(
        '--vcfFile',
        help=('The VCF file. If omitted, bcftools will be used to make a VCF '
              'file from the BAM file.'))

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '--id',
        help=('The id to use in the consensus sequence in the output FASTA. '
              'If not given, the reference sequence id will be used.'))

    group.add_argument(
        '--idLambda', metavar='LAMBDA-FUNCTION',
        help=('A one-argument function taking and returning a read id. '
              'This can be used to set the id of the reference sequence based '
              'on the id of the reference sequence (the function will be '
              'called with the id of the reference sequence). E.g., '
              '--idLambda "lambda id: id.split(\'_\')[0]" or '
              '--idLambda "lambda id: id[:10] + \'-consensus\'".'))

    parser.add_argument(
        '--sample',
        help=('The name of the sample (from the @RG SM tag in the original '
              'alignment BAM file) for which a consensus should be made. '
              'If not given, the first sample name (from the #CHROM header) '
              'in the VCF file will be used.'))

    parser.add_argument(
        '--dryRun', default=False, action='store_true',
        help='Do not run commands, just print what would be done.')

    parser.add_argument(
        '--maskLowCoverage', default=0, type=int,
        help=('Put an N into sites where the coverage is below the specified '
              'cutoff. If you specify a negative numer, masking will be '
              'turned off. Requires --bam.'))

    parser.add_argument(
        '--log', default=False, action='store_true',
        help=('Show a log of commands that were (or would be, if --dryRun is '
              'used) executed.'))

    parser.add_argument(
        '--noClean', default=True, action='store_false', dest='clean',
        help=('Do not remove intermediate files or the temporary directory.'))

    parser.add_argument(
        '--callHaplotypesGATK', default=False, action='store_true',
        help=('Use GATK to call haplotypes. See '
              'https://gatk.broadinstitute.org for details on GATK.'))

    parser.add_argument(
        '--picardJar',
        help=('The path to the Picard jar file. See '
              'https://github.com/broadinstitute/picard for details on '
              'Picard.'))

    parser.add_argument(
        '--ivar', default=False, action='store_true',
        help='If given, ivar will be used to call the consensus.')

    parser.add_argument(
        '--ivarFrequencyThreshold', type=float,
        help=(f'The frequency threshold used by ivar when calling the '
              f'consensus. If the frequency of the most-common nucleotide at '
              f'a site meets this threshold, the nucleotide will be called. '
              f'Otherwise, an ambiguous nucleotide code will be produced, '
              f'based on the smallest set of most-frequent nucleotides whose '
              f'summed frequencies meet the threshold. See {IVAR_DOCS} for '
              f'more information. If not given, '
              f'{IVAR_FREQUENCY_THRESHOLD_DEFAULT} is used. Can only be used '
              f'if --ivar is also specified.'))

    parser.add_argument(
        '--ivarBedFile',
        help=('If ivar should trim primers, a bed file of the primer '
              'positions.'))

    args = parser.parse_args()

    if not (args.bam or args.vcfFile):
        print('At least one of --bam or --vcfFile must be given.',
              file=sys.stderr)
        sys.exit(1)

    if args.maskLowCoverage and not args.bam:
        print('If --maskLowCoverage is used, --bam must be too.',
              file=sys.stderr)
        sys.exit(1)

    if args.ivar and not args.bam:
        print('If --ivar is used, --bam must be too.', file=sys.stderr)
        sys.exit(1)

    if args.ivarFrequencyThreshold is not None and not args.ivar:
        print('If --ivarFrequencyThreshold is used, --ivar must be too.',
              file=sys.stderr)
        sys.exit(1)

    if args.ivar and args.ivarFrequencyThreshold is None:
        args.ivarFrequencyThreshold = IVAR_FREQUENCY_THRESHOLD_DEFAULT

    e = Executor(args.dryRun)

    tempdir = mkdtemp(prefix='consensus-')

    if args.vcfFile:
        vcfFile = args.vcfFile
    else:
        # No VCF file provided, so make one.
        vcfFile = join(tempdir, 'vcf.gz')
        if args.callHaplotypesGATK:
            e.execute("samtools index '%s'" % args.bam)
            if args.picardJar:
                picardJar = args.picardJar
            else:
                try:
                    picardJar = os.environ['PICARD_JAR']
                except KeyError:
                    print('If you use --callHaplotypesGATK, you must give a '
                          'Picard JAR file with --picardJar or else set '
                          'PICARD_JAR in your environment.', file=sys.stderr)
                    sys.exit(1)

            indexFile = args.reference + '.fai'
            if os.path.exists(indexFile):
                removeIndex = False
            else:
                removeIndex = True
                e.execute("samtools faidx '%s'" % args.reference)

            if args.reference.lower().endswith('.fasta'):
                dictFile = args.reference[:-len('.fasta')] + '.dict'
            else:
                dictFile = args.reference + '.dict'

            if os.path.exists(dictFile):
                removeDict = False
            else:
                removeDict = True
                e.execute(
                    "java -jar '%s' CreateSequenceDictionary R='%s' O='%s'"
                    % (picardJar, args.reference, dictFile))

            e.execute(
                'gatk --java-options -Xmx4g HaplotypeCaller '
                "--reference '%s' "
                "--input '%s' "
                "--output '%s' "
                "--sample-ploidy 1 "
                '-ERC GVCF' %
                (args.reference, args.bam, vcfFile))

            if removeIndex:
                e.execute("rm '%s'" % indexFile)

            if removeDict:
                e.execute("rm '%s'" % dictFile)
        else:
            e.execute("bcftools mpileup --max-depth 5000 -Ou -f '%s' '%s' | "
                      "bcftools call --ploidy 1 -mv -Oz -o '%s'" %
                      (args.reference, args.bam, vcfFile))

            e.execute("bcftools index '%s'" % vcfFile)

    if args.maskLowCoverage >= 0:
        # Make a BED file.
        bedFile = join(tempdir, 'mask.bed')
        # The doubled-% below are so that Python doesn't try to fill in the
        # values and instead just generates a single % that awk sees.
        e.execute(
            "samtools depth -a '%s' | "
            "awk '$3 < %d {printf \"%%s\\t%%d\\t%%d\\n\", "
            "$1, $2 - 1, $2}' > '%s'" %
            (args.bam, args.maskLowCoverage, bedFile))
        maskArg = '--mask ' + bedFile
    else:
        maskArg = ''

    if args.sample:
        sample = args.sample
    else:
        result = e.execute(
            "gunzip -c '%s' | egrep -m 1 '^#CHROM' | cut -f10" % vcfFile)
        sample = result.stdout.strip()

    consensusFile = join(tempdir, 'consensus.fasta')

    if args.ivar:
        if args.ivarBedFile:
            tempBamFile = join(tempdir, basename(args.bam) + '-trimmed')
            result = e.execute(
                "ivar trim -i %r -b %r -p %r -e" % (
                    args.bam, args.ivarBedFile, tempBamFile))
            ivarTempBamFile = tempBamFile + '.bam'
            sortedIvarTempBamFile = tempBamFile + '-trimmed-sorted.bam'
            result = e.execute(
                "samtools sort %r -o %r" % (
                    ivarTempBamFile, sortedIvarTempBamFile))
            bamFile = sortedIvarTempBamFile
        else:
            bamFile = args.bam

        ivarConsensusFile = join(tempdir, 'temporary-consensus')
        result = e.execute(
            "samtools mpileup -A -Q 0 %r | "
            "ivar consensus -p %r -q 20 -t %r -m %r" % (
                bamFile, ivarConsensusFile, args.ivarFrequencyThreshold,
                args.maskLowCoverage))

        result = e.execute(
            "mv %s %s" % (ivarConsensusFile + '.fa', consensusFile))

    else:
        result = e.execute(
            "bcftools consensus --sample '%s' --iupac-codes %s --fasta-ref "
            "'%s' '%s' > '%s'" %
            (sample, maskArg, args.reference, vcfFile, consensusFile))

    consensus = list(FastaReads(consensusFile))[0]
    if args.id is not None:
        consensus.id = args.id
    elif args.idLambda is not None:
        idLambda = eval(args.idLambda)
        consensus.id = idLambda(consensus.id)

    print(consensus.toString('fasta'), end='')

    if result.stderr:
        print(result.stderr, end='', file=sys.stderr)

    if args.dryRun or args.log:
        print('\n'.join(e.log), file=sys.stderr)

    if tempdir:
        if args.clean:
            e.execute("rm -r '%s'" % tempdir)
        else:
            print('Temporary directory %r.' % tempdir, file=sys.stderr)


if __name__ == '__main__':
    main()

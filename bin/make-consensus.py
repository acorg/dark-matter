#!/usr/bin/env python

import sys
import argparse
from tempfile import mkdtemp
from os.path import join

from dark.fasta import FastaReads
from dark.process import Executor


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
              'Required if --maskNoCoverage is used. If no BAM file is '
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
        '--maskNoCoverage', default=False, action='store_true',
        help='Put an N into sites where there is no coverage. Requires --bam.')

    parser.add_argument(
        '--log', default=False, action='store_true',
        help=('Show a log of commands that were (or would be, if --dryRun is '
              'used) executed.'))

    parser.add_argument(
        '--noClean', default=True, action='store_false', dest='clean',
        help=('Do not remove intermediate files or the temporary directory.'))

    args = parser.parse_args()

    if not (args.bam or args.vcfFile):
        print('At least one of --bam or --vcfFile must be given.',
              file=sys.stderr)
        sys.exit(0)

    if args.maskNoCoverage and not args.bam:
        print('If --maskNoCoverage is used, --bam must be too.',
              file=sys.stderr)
        sys.exit(0)

    e = Executor(args.dryRun)

    tempdir = mkdtemp(prefix='consensus-')

    if args.vcfFile:
        vcfFile = args.vcfFile
    else:
        # No VCF file provided, so make one.
        vcfFile = join(tempdir, 'vcf.gz')
        e.execute("bcftools mpileup --max-depth 5000 -Ou -f '%s' '%s' | "
                  "bcftools call -mv -Oz -o '%s'" %
                  (args.reference, args.bam, vcfFile))

        e.execute("bcftools index '%s'" % vcfFile)

    if args.maskNoCoverage:
        # Make a BED file.
        bedFile = join(tempdir, 'mask.bed')
        # The doubled-% below are so that Python doesn't try to fill in the
        # values and instead just generates a single % that awk sees.
        e.execute(
            "samtools depth -a '%s' | "
            "awk '$3 == 0 {printf \"%%s\\t%%d\\t%%d\\n\", "
            "$1, $2 - 1, $2}' > '%s'" %
            (args.bam, bedFile))
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

#!/usr/bin/env python

import sys
import argparse
from tempfile import mkdtemp
from os.path import join

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
        help=('The BAM file from which the consensus should be made. Only '
              'used if a VCF file is not provided.'))

    parser.add_argument(
        '--sample',
        help=('The name of the sample (from the @RG SM tag in the original '
              'alignment BAM file for which a consensus should be made. '
              'If not given, the first sample name in the VCF file will be '
              'used.'))

    parser.add_argument(
        '--vcfFile',
        help=('The VCF file. If ommitted, bcftools will be used to make a VCF '
              'file.'))

    parser.add_argument(
        '--dryRun', default=False, action='store_true',
        help='Do not run commands, just print what would be done.')

    parser.add_argument(
        '--log', default=False, action='store_true',
        help=('Show a log of commands that were (or would be, if --dryRun is '
              'used) executed.'))

    parser.add_argument(
        '--noClean', default=True, action='store_false', dest='clean',
        help=('Do not remove intermediate files or the temporary directory.'))

    args = parser.parse_args()

    if args.bam and args.vcfFile:
        print('Only of --bam or --vcfFile can be given.', file=sys.stderr)
        sys.exit(0)

    if not (args.bam or args.vcfFile):
        print('One of --bam or --vcfFile must be given.', file=sys.stderr)
        sys.exit(0)

    e = Executor(args.dryRun)

    if args.bam:
        # Make a VCF file.
        tempdir = mkdtemp(prefix='consensus-')
        vcfFile = join(tempdir, 'vcf.gz')
        e.execute("bcftools mpileup -Ou -f '%s' '%s' | "
                  "bcftools call -mv -Oz -o '%s'" %
                  (args.reference, args.bam, vcfFile))

        e.execute("bcftools index '%s'" % vcfFile)
    else:
        tempdir = None
        vcfFile = args.vcfFile

    if args.sample:
        sample = args.sample
    else:
        result = e.execute("gunzip -c '%s' | egrep '^#CHROM' | cut -f10" %
                           vcfFile)

        sample = result.stdout.split('\n')[0]

    result = e.execute(
        "bcftools consensus --sample '%s' --fasta-ref '%s' '%s'" %
        (sample, args.reference, vcfFile))

    print(result.stdout, end='')
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

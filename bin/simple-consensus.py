#!/usr/bin/env python

import sys
import argparse

from dark.consensus import consensusFromBAM
from dark.fasta import FastaReads
from dark.reads import DNARead


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Make a consensus sequence.')

    parser.add_argument(
        '--bamFilename', required=True,
        help='The BAM file from which the consensus should be called.')

    parser.add_argument(
        '--referenceId',
        help=('The id of the reference sequence in the --bam file. '
              'If not given, the reference id (as given by --referenceId or '
              'taken from the first FASTA record in the --reference file) '
              'will be used. Caution should be taken when using this option '
              'as it (deliberately) allows the names of the FASTA reference '
              'sequence and the reference referred to in the BAM file to '
              'differ. This can be convenient if a sequence with an '
              'identical name (to that of the BAM file) is not present in the '
              'reference FASTA file. But care must be taken because this '
              'allows you to accidentally make a consensus using a reference '
              'sequence that is not the one that was used to make the BAM '
              'file. The length of the reference sequence (if a reference '
              'sequence is given) must match the length of the '
              'reference mentioned in the BAM file. If no reference or BAM '
              'id is given, a consensus will only be called if the BAM file '
              'refers to a single reference.'))

    parser.add_argument(
        '--referenceFasta',
        help='The reference FASTA file.')

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '--consensusId',
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
        '--minCoverage', default=0, type=int,
        help=('The minimum number of reads that must cover a site for a '
              'consensus base to be called. If zero reads cover a site, the '
              '--noCoverage value is used or if the number is greater than '
              'zero but less than minCoverage, the --lowCoverage value is '
              'used.'))

    parser.add_argument(
        '--lowCoverage', default='reference',
        help=('What to do when some reads cover a site, but fewer than '
              'the --minCoverage value. Either set this to "reference" (to '
              'use the base from the reference sequence)or a single character '
              '(e.g., "N").'))

    parser.add_argument(
        '--noCoverage', default='reference',
        help=('What to do when no reads cover a site, Either set this to '
              '"reference" (to use the base from the reference sequence), or '
              'a single character (e.g., "N").'))

    parser.add_argument(
        '--threshold', type=float, default=0.6,
        help=('The frequency threshold when calling the consensus. If the '
              'frequency of the most-common nucleotide at a site meets this '
              'threshold, that nucleotide will be called. Otherwise, an '
              'ambiguous nucleotide code will be produced, based on the '
              'smallest set of most-frequent nucleotides whose summed '
              'frequencies meet the threshold. If the frequency of the '
              'nucleotide that causes the threshold to be reached is the '
              'same as that of other nucleotides, all such nucleotides will '
              'be included in the ambiguous code.'))

    parser.add_argument(
        '--ignoreQuality', action='store_true',
        help=('Ignore FASTQ quality scores. All bases of all reads will be '
              'considered equally reliable.'))

    parser.add_argument(
        '--includeSoftClipped', action='store_true',
        help=('Include information from read bases that were marked as '
              'soft-clipped by the algorithm that made the BAM file. See '
              'https://samtools.github.io/hts-specs/SAMv1.pdf if this makes '
              'no sense to you.'))

    parser.add_argument(
        '--compareWithPileupFile',
        help=('Make a consensus using the pysam pileup function and compare '
              'it to the one made by usign pysam fetch. Write the summary to '
              'this file.'))

    parser.add_argument(
        '--deletionSymbol', default='-',
        help=('When a deletion is detected, put this symbol into the '
              'consensus sequence. Use the empty string to have deletions '
              'omitted.'))

    parser.add_argument(
        '--strategy', default='fetch', choices=('fetch',),
        help='The consensus-making strategy to use.')

    parser.add_argument(
        '--deletionThreshold', default=0.5, type=float,
        help=('If some reads have a deletion at a site and some do not, call '
              'the site as a deletion if the fraction of reads with the '
              'deletion is at least this value.'))

    parser.add_argument(
        '--insertionCountThreshold', default=5, type=int,
        help=('The number of reads that must have an insertion at an offset '
              'in order for the insertion to be called in the consensus.'))

    parser.add_argument(
        '--progress', action='store_true',
        help='Show a progress bar.')

    args = parser.parse_args()

    if args.referenceFasta is None:
        if args.lowCoverage == 'reference':
            print('The --lowCoverage option is set to "reference" '
                  'but no reference sequence has been supplied',
                  file=sys.stderr)
            sys.exit(1)
        if args.noCoverage == 'reference':
            print('The --noCoverage option is set to "reference" '
                  'but no reference sequence has been supplied',
                  file=sys.stderr)
            sys.exit(1)
        reference = None
    else:
        if args.referenceId is None:
            # Use the first sequence as the reference.
            reference = list(FastaReads(args.referenceFasta))[0]
            args.referenceId = reference.id
        else:
            for read in FastaReads(args.referenceFasta):
                if read.id == args.referenceId:
                    reference = read
                    break
            else:
                print(f'Could not find a sequence with id '
                      f'{args.referenceId!r} in {args.referenceFasta!r}',
                      file=sys.stderr)
                sys.exit(1)

    if args.consensusId is not None:
        consensusId = args.consensusId
    elif args.idLambda is not None:
        if reference is None:
            print('You used --idLambda, but no reference was given (so '
                  'there is no id to apply the lambda function to).',
                  file=sys.stderr)
            sys.exit(1)
        idLambda = eval(args.idLambda)
        consensusId = idLambda(reference.id)
    else:
        consensusId = ((args.referenceId or 'consensus') if reference is None
                       else reference.id)

    consensus = consensusFromBAM(
        args.bamFilename, referenceId=args.referenceId, reference=reference,
        threshold=args.threshold, minCoverage=args.minCoverage,
        lowCoverage=args.lowCoverage, noCoverage=args.noCoverage,
        deletionSymbol=args.deletionSymbol,
        deletionThreshold=args.deletionThreshold,
        insertionCountThreshold=args.insertionCountThreshold,
        ignoreQuality=args.ignoreQuality,
        includeSoftClipped=args.includeSoftClipped, strategy=args.strategy,
        compareWithPileupFile=args.compareWithPileupFile,
        progress=args.progress)

    print(DNARead(consensusId, consensus).toString('fasta'), end='')


if __name__ == '__main__':
    main()

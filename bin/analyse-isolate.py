#!/usr/bin/env python

from Bio.Data.CodonTable import ambiguous_dna_by_id

from dark.fasta import FastaReads
from dark.reads import Reads
from dark.sam import SAMFilter

from mvlib.minorVariants import MinorVariantInfo
from mvlib.functions import isMinorVariantPosition

from sars2seq.alignment import SARS2Alignment, offsetInfoMultipleGenomes
from sars2seq.features import Features


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Analyse isolate sequencing data.'))

    parser.add_argument(
        '--alignmentFile',
        help=('The name of the alignment file. Either a bam file or a json '
              'file that has been generated using MinorVariantInfo.'))

    parser.add_argument(
        '--clinicalSequence',
        help=('The fasta file of the sequence of the clinical isolate that '
              'the isolate was made from.'))

    parser.add_argument(
        '--consensusSequence',
        help=('The fasta file of the consensus sequence from the isolate '
              'sequencing run.'))

    parser.add_argument(
        '--sampleName',
        help='The name of the sample.')

    parser.add_argument(
        '--minorVariantFrequencyCutoff', default=0.05, type=float,
        help='The cut-off for minor variant frequency.')

    parser.add_argument(
        '--minorVariantCoverageCutoff', default=5, type=int,
        help=('The minimum number of reads required to evaluate a sequence '
              'for minority variants.'))

    args = parser.parse_args()

    codonTable = ambiguous_dna_by_id[1].forward_table

    # Make the genome objects for the clinical sequence and the consensus
    # sequence of the isolate sequencing run.
    clinicalSeq = list(FastaReads(args.clinicalSequence))[0]
    consensusSeq = list(FastaReads(args.consensusSequence))[0]
    features = Features()
    clinicalAlignment = SARS2Alignment(clinicalSeq, features)
    consensusAlignment = SARS2Alignment(consensusSeq, features)

    # Get the minor variant information
    if args.alignmentFile.endswith('.json'):
        mvInfo = MinorVariantInfo(jsonFile=args.alignmentFile)
    else:
        mvInfo = MinorVariantInfo(
            bamFile=args.alignmentFile,
            minBaseQuality=35, minMappingQuality=30)
        referenceLengths = SAMFilter(
            args.alignmentFile,
            filterRead=Reads().filter().filterRead).referenceLengths()
        assert 'hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125' in referenceLengths, (
               'This looks like an incorrect reference')

    coverageDepthSum = 0
    coverageDepthMin = 10000
    coverageDepthMax = 0

    print('Isolate analysis for sample %s' % args.sampleName)
    print('Sequence ID of clinical sample: %s' % clinicalSeq.id)

    for position in range(len(consensusSeq)):
        # 'position' here will be relative to the reference sequence (has to
        # be EPI_ISL_402125/NC_045512.2) that the reads from the isolate were
        # aligned against.

        # Count the coverage
        readCountAtPos = sum(mvInfo.countsPerBase[position].values())
        coverageDepthSum += readCountAtPos
        if coverageDepthMin > readCountAtPos:
            coverageDepthMin = readCountAtPos
        if coverageDepthMax < readCountAtPos:
            coverageDepthMax = readCountAtPos

        # Get just one feature at this position
        feats = features.getFeatureNames(position, includeUntranslated=True)

        feature = sorted(list(feats))[0] if feats else None

        offsetInfo = offsetInfoMultipleGenomes(
            [clinicalAlignment, consensusAlignment], position,
            featureName=feature, relativeToFeature=False,
            includeUntranslated=True)

        clinInfo = offsetInfo['genomes'][0]
        clinBase = clinInfo['codon'][clinInfo['frame']]
        consInfo = offsetInfo['genomes'][1]
        consBase = consInfo['codon'][consInfo['frame']]

        # Check if there are any differences between the consensus from
        # the isolate and the consensus from the clinical sample.
        if clinBase != consBase:
            furin = ''
            if clinInfo['aa'] == consInfo['aa']:
                synInfo = 'synonymous'
            else:
                synInfo = 'non-synonymous (%s, %s%d%s)' % (
                    offsetInfo['featureName'], clinInfo['aa'],
                    offsetInfo['reference']['aaOffset'] + 1, consInfo['aa'])
                if (offsetInfo['featureName'] == 'surface glycoprotein' and
                        670 < offsetInfo['reference']['aaOffset'] < 690):
                    furin = ' FURIN SITE!'
            if clinBase != '-' and 70 < position < 29714:
                print('\tChange from clinical sample sequence at position %d '
                      '(nt change: %s>%s), %s%s' % (
                          position + 1, clinBase, consBase, synInfo, furin))

        # Check if there are any minority variants at this position.
        if isMinorVariantPosition(mvInfo.countsPerBase[position],
                                  args.minorVariantCoverageCutoff,
                                  args.minorVariantFrequencyCutoff):

            clinCod = clinInfo['codon']
            consCod = consInfo['codon']

            print('\t\tMinority variant at position: %d (%s), '
                  'coverage: %d reads' % (
                      position + 1, offsetInfo['featureName'],
                      sum(mvInfo.countsPerBase[position].values())))
            for base in mvInfo.countsPerBase[position]:
                if (mvInfo.countsPerBase[position][base] >
                        args.minorVariantCoverageCutoff):
                    mvCod = (
                        consInfo['codon'][:consInfo['frame']] +
                        base + consInfo['codon'][consInfo['frame'] + 1:])
                    if clinBase == base:
                        baseInfo = 'consensus base'
                    else:
                        # Get the amino acid
                        if 'N' in clinCod:
                            clinCodT = 'X'
                        else:
                            clinCodT = codonTable.get(clinCod)
                        if 'N' in mvCod:
                            mvCodT = 'X'
                        else:
                            mvCodT = codonTable.get(mvCod)
                        # Compare the amino acids
                        if clinCodT == mvCodT:
                            baseInfo = 'synonymous'
                        else:
                            baseInfo = '%s%d%s, %s' % (
                                clinCodT,
                                offsetInfo['reference']['aaOffset'] + 1,
                                mvCodT, offsetInfo['featureName'])

                    print('\t\t\tBase: %s, count: %d (%.2f): %s' % (
                        base, mvInfo.countsPerBase[position][base],
                        mvInfo.countsPerBase[position][base] / sum(
                            mvInfo.countsPerBase[position].values()),
                        baseInfo))

    print('Genome coverage: %.2f, mean coverage depth: %.2f (min: %d, '
          'max: %d)\n' % (
              (len(consensusSeq) -
               consensusSeq.sequence.count('N')) / len(consensusSeq),
              coverageDepthSum / len(consensusSeq),
              coverageDepthMin, coverageDepthMax))

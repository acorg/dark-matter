# from time import time
from collections import defaultdict

from dark.dna import leastAmbiguousFromCounts
from dark.sam import samfile


def consensusFromBAM(bamFilename, reference, strategy='majority',
                     threshold=0.8, minCoverage=1, lowCoverage='reference',
                     noCoverage='reference'):
    """
    Build a consensus sequence from a BAM file.

    @param bamFilename: the BAM file.
    @param reference: A C{Read} instance giving the reference sequence.
    @param strategy: A C{str} strategy, one of 'majority'.
    @param threshold: A C{float} threshold. This fraction, at least, of the
        most-common nucleotides at a site are used to determine the consensus
        nucleotide (or ambiguous symbol if more than one nucleotide is
        required to achieve this threshold). If there is a tie in nucleotide
        counts at a site that causes the threshold to be met, all nucleotides
        of equeal frequncy will be included in the ambiguous symbol for that
        site. This is perhaps better explained with an example. See
        https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
        and the corresponding testGeneiousExamplesTie test in test/test_dna.py
    @param minCoverage: An C{int} minimum number of reads that must cover a
        site for a threshold consensus base to be called. If zero reads
        cover a site, the C{noCoverage} value is used or if the number is
        greater than zero but less then C{minCoverage}, the C{lowCoverage}
        value is used.
    @param lowCoverage: A C{str} indicating what to do when some reads cover a
        site, but fewer than C{minCoverage}. Either 'reference' or a single
        character (e.g., 'N').
    @param noCoverage: A C{str} indicating what to do when no reads cover a
        reference base. Either 'reference' or a single character (e.g., 'N').
    @return: A C{str} consensus sequence.
    """
    with samfile(bamFilename) as fp:
        if strategy == 'majority':
            return _majorityConsensus(fp, reference, threshold, minCoverage,
                                      lowCoverage, noCoverage)
        else:
            raise ValueError(f'Unknown consensus strategy {strategy!r}.')


def _majorityConsensus(bam, reference, threshold, minCoverage, lowCoverage,
                       noCoverage):
    """
    Compute a majority consensus.

    @param bam: An open BAM file.
    @param reference: A C{dark.reads.Read} instance giving the reference
        sequence.
    @param threshold: A C{float} threshold. If the majority ...
    @param minCoverage: An C{int} minimum number of reads that must cover a
        site for a threshold consensus base to be called. If zero reads
        cover a site, the C{noCoverage} value is used or if the number is
        greater than zero but less then C{minCoverage}, the C{lowCoverage}
        value is used.
    @param lowCoverage: A C{str} indicating what to do when some reads cover a
        site, but fewer than C{minCoverage}. Either 'reference' or a single
        character (e.g., 'N').
    @param noCoverage: A C{str} indicating what to do when no reads cover a
        reference base. Either 'reference' or a single character (e.g., 'N').
    @return: A C{str} consensus sequence.
    """
    result = list(reference.sequence if noCoverage == 'reference' else
                  noCoverage * len(reference))
    lowCoverage = (reference.sequence if lowCoverage == 'reference' else
                   lowCoverage * len(reference))

    for column in bam.pileup(reference=reference.id):
        site = column.reference_pos
        bases = defaultdict(int)
        readCount = 0
        for read in column.pileups:
            readCount += 1
            base = read.alignment.query_sequence[read.query_position]
            bases[base] += 1
            # quality = read.alignment.query_qualities[read.query_position]

        result[site] = (lowCoverage[site] if readCount < minCoverage else
                        leastAmbiguousFromCounts(bases, threshold))

    return ''.join(result)

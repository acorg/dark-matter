import sys
# from time import time
from collections import defaultdict

from dark.dna import leastAmbiguousFromCounts
from dark.sam import (
    samfile, samReferences, UnequalReferenceLengthError,
    UnknownReference, UnspecifiedReference)


def consensusFromBAM(bamFilename, referenceId=None, reference=None,
                     threshold=0.8, minCoverage=1, lowCoverage='reference',
                     noCoverage='reference', ignoreQuality=False,
                     strategy='majority', logfp=None):
    """
    Build a consensus sequence from a BAM file.

    @param bamFilename: the BAM file.
    @param referenceId: A C{str} reference name indicating which reference to
        find reads against in the BAM file. If C{None} and there is only one
        reference in the BAM file, that one will be used, else a ValueError
        is raised due to not knowing which reference to use.
    @param reference: A C{Read} instance giving the reference sequence, or
        C{None} (in which case C{lowCoverage} and C{noCoverage} may not be
        'reference'.
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
        site for a consensus base to be called. If zero reads cover a site, the
        C{noCoverage} value is used or if the number is greater than zero but
        less than C{minCoverage}, the C{lowCoverage} value is used.
    @param lowCoverage: A C{str} indicating what to do when some reads cover a
        site, but fewer than C{minCoverage}. Either 'reference' or a single
        character (e.g., 'N').
    @param noCoverage: A C{str} indicating what to do when no reads cover a
        reference base. Either 'reference' or a single character (e.g., 'N').
    @param ignoreQuality: If C{True}, ignore quality scores.
    @param strategy: A C{str} consensus-making strategy (cuurently must be
        'majority').
    @raise UnspecifiedReference: If no id is provided to indicate which BAM
        file reference to call a consensus for.
    @raise UnknownReference: If a requested reference id is unknown.
    @raise UnequalReferenceLengthError: If the passed reference does not have a
        length identical to the length mentioned in the BAM file.
    @return: A C{str} consensus sequence.
    """
    with samfile(bamFilename) as bam:
        bamReferences = set(samReferences(bam))
        if referenceId is None:
            if reference:
                referenceId = reference.id
            elif len(bamReferences) == 1:
                referenceId = tuple(bamReferences)[0]

        if referenceId is None:
            raise UnspecifiedReference(
                f'BAM file {str(bamFilename)!r} mentions '
                f'{len(bamReferences)} references '
                f'({", ".join(sorted(bamReferences))}) but you have not '
                f'passed a referenceId argument or a reference sequence to '
                f'indicate which one to use.')

        tid = bam.get_tid(referenceId)

        if tid == -1 or referenceId not in bamReferences:
            raise UnknownReference(
                f'BAM file {str(bamFilename)!r} does not mention a '
                f'reference with id {referenceId!r}.')

        referenceLength = bam.lengths[tid]

        if reference and len(reference) != referenceLength:
            raise UnequalReferenceLengthError(
                f'Reference with id {reference.id!r} has length '
                f'{len(reference)}, which does not match the length of '
                f'reference {referenceId!r} ({referenceLength}) in BAM file '
                f'{str(bamFilename)!r}.')

        if strategy == 'majority':
            return _majorityConsensus(
                bam, referenceId, reference, referenceLength, threshold,
                minCoverage, lowCoverage, noCoverage, ignoreQuality, logfp)
        else:
            raise ValueError(f'Unknown consensus strategy {strategy!r}.')


def _majorityConsensus(bam, referenceId, reference, referenceLength, threshold,
                       minCoverage, lowCoverage, noCoverage, ignoreQuality,
                       logfp):
    """Compute a majority consensus.

    @param bam: An open BAM file.
    @param referenceId: A C{str} reference name indicating which reference to
        find reads against in the BAM file. If C{None} and there is only one
        reference in the BAM file, that one will be used, else a ValueError
        is raised due to not knowing which reference to use.
    @param reference: A C{Read} instance giving the reference sequence, or
        C{None} (in which case neither C{lowCoverage} nor C{noCoverage} may be
        'reference').
    @param referenceLength: The C{int} length of the reference (note that we
        may not have the reference sequence but we can still get its length
        from the BAM file header).
    @param threshold: A C{float} threshold used when calling the consensus.
        If the frequency of the most-common nucleotide at a site meets this
        threshold, that nucleotide will be called. Otherwise, an ambiguous
        nucleotide code will be produced, based on the smallest set of
        most-frequent nucleotides whose summed frequencies meet the
        threshold. If the frequency of the nucleotide that causes the
        threshold to be reached is the same as that of other nucleotides,
        all such nucleotides will be included in the ambiguous code.
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
    @param ignoreQuality: If C{True}, ignore quality scores.
    @param logfp: If not C{None}, an open file pointer for writing information
        to.
    @return: A C{str} consensus sequence.

    """
    result = list(reference.sequence if noCoverage == 'reference' else
                  noCoverage * referenceLength)
    lowCoverage = (reference.sequence if lowCoverage == 'reference' else
                   lowCoverage * referenceLength)

    deletedSites = defaultdict(int)

    for column in bam.pileup(reference=referenceId):
        site = column.reference_pos
        # print(f'COL {column.reference_pos}', file=sys.stderr)
        bases = defaultdict(int)
        readCount = 0
        for readCount, read in enumerate(column.pileups, start=1):
            # print(f'  qp {read.query_position!r}')
            if read.is_del:
                # print(f'DEL at site {site}', read, file=sys.stderr)
                deletedSites[site] += 1
                continue
            elif read.is_refskip:
                # The read skips a site in the reference. I.e., there is
                # a deletion in this read.
                # print('REF_SKIP:', read, file=sys.stderr)
                raise NotImplementedError()
            else:
                print('NORM:', read, file=sys.stderr)
                queryPosition = read.query_position
                base = read.alignment.query_sequence[queryPosition]

                bases[base] += (1 if ignoreQuality else
                                read.alignment.query_qualities[queryPosition])

        result[site] = (lowCoverage[site] if readCount < minCoverage else
                        leastAmbiguousFromCounts(bases, threshold))

    newResult = []
    for site, base in enumerate(result):
        if site not in deletedSites:
            newResult.append(base)

    if logfp:
        for site, count in deletedSites.items():
            print(f'{site}: {count} deletions', file=logfp)

    return ''.join(newResult)

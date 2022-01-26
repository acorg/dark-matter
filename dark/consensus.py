import sys
from time import time
from collections import defaultdict
from contextlib import contextmanager
import progressbar

from dark.dna import leastAmbiguousFromCounts
from dark.sam import (
    samfile, samReferences, UnequalReferenceLengthError,
    UnknownReference, UnspecifiedReference)

DEBUG = False
# DEBUG = True


def debug(*msg):
    print(*msg, file=sys.stderr)


class Insertions:
    def __init__(self, n):
        self.n = n
        self.sets = [set() for _ in range(n)]

    def add(self, i):
        self.sets[i % self.n].add(i)

    def countLE(self, value):
        start = time()
        count = 0
        for s in self.sets:
            for i in s:
                if i <= value:
                    count += 1
        stop = time()
        total = sum([len(s) for s in self.sets])
        if DEBUG:
            debug(f'Counting in {total:5d} insert indices took '
                  f'{stop - start:.6f}  Result={count:5d}')

        return count


class Bases:
    __slots__ = ('count', 'counts')

    def __init__(self):
        self.count = 0
        # self.counts = defaultdict(int)
        self.counts = dict.fromkeys('ACGT', 0)

    def __str__(self):
        return f'<Bases count={self.count}, bases={self.counts}'

    def __repr__(self):
        return str(self)

    def add(self, base, quality):
        if base != 'N':
            self.count += 1
            self.counts[base] += quality


@contextmanager
def maybeProgressBar(show, nReads):
    """
    A context manager to maybe show a progress bar.

    @param show: If C{True}, yield a progress bar, else yield C{None}.
    @param nReads: The C{int} number of mapped reads that will be processed.
    """
    if show:
        with progressbar.ProgressBar(max_value=nReads) as bar:
            yield bar
    else:
        yield None


def consensusFromBAM(bamFilename, referenceId=None, reference=None,
                     threshold=0.8, minCoverage=1, lowCoverage='reference',
                     noCoverage='reference', ignoreQuality=False,
                     strategy='majority', logfp=None, progress=False):
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
    @param progress: If C{True}, display a progress bar on standard error.
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
                f'reference with id {referenceId!r}. Known references are: '
                f'{", ".join(sorted(bamReferences))}.')

        referenceLength = bam.lengths[tid]

        if reference and len(reference) != referenceLength:
            raise UnequalReferenceLengthError(
                f'Reference with id {reference.id!r} has length '
                f'{len(reference)}, which does not match the length of '
                f'reference {referenceId!r} ({referenceLength}) in BAM file '
                f'{str(bamFilename)!r}.')

        if strategy == 'majority':
            return _fetchConsensus(
                bam, referenceId, reference, referenceLength, threshold,
                minCoverage, lowCoverage, noCoverage, ignoreQuality, progress,
                logfp)
        else:
            raise ValueError(f'Unknown consensus strategy {strategy!r}.')


def _fetchConsensus(bam, referenceId, reference, referenceLength, threshold,
                    minCoverage, lowCoverage, noCoverage, ignoreQuality,
                    progress, logfp):
    """Compute a majority consensus using fetch.

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
    correspondences = defaultdict(lambda: Bases())
    deletions = set()
    insertions = defaultdict(lambda: Bases())
    insertionSet = Insertions(1)

    nReads = bam.count(contig=referenceId)
    with maybeProgressBar(progress, nReads) as bar:
        for readCount, read in enumerate(bam.fetch(contig=referenceId),
                                         start=1):

            addPairsInfo(
                read.get_aligned_pairs(), read.query_sequence,
                ([1] * len(read.query_sequence) if ignoreQuality else
                 read.query_qualities),
                referenceLength, correspondences, deletions, insertions,
                insertionSet)

            if DEBUG:
                debug(f'read id  : {read.query_name}')
                debug('query    :', read.query_sequence)
                debug(f'cigar    : {read.cigarstring}')
                debug(f'match    : {read.reference_start}')
                debug(f'Pairs    : {read.get_aligned_pairs()}')
                debug(f'  {correspondences=}')
                debug(f'  {insertions=}')
                debug(f'  {deletions=}')

            if bar:
                bar.update(readCount)

    orig = list(reference.sequence if noCoverage == 'reference' else
                noCoverage * referenceLength)
    lowCoverageStr = (reference.sequence if lowCoverage == 'reference' else
                      lowCoverage * referenceLength)

    minCorrespondence = min(correspondences, default=0)
    maxCorrespondence = max(correspondences, default=0)

    prefix = [None] * (abs(minCorrespondence) if minCorrespondence < 0 else 0)
    suffix = [None] * (maxCorrespondence - referenceLength + 1
                       if maxCorrespondence >= referenceLength else 0)

    result = list(orig)
    for offset, bases in correspondences.items():
        if offset < 0:
            array = prefix
            offset += len(prefix)
        elif offset >= referenceLength:
            array = suffix
            offset = offset - referenceLength
        else:
            array = result

        if bases.count < minCoverage:
            array[offset] = lowCoverageStr[offset]
        else:
            array[offset] = leastAmbiguousFromCounts(bases.counts, threshold)

    # Do deletions.
    resultWithDeletions = []
    for offset, base in enumerate(result):
        if offset not in deletions:
            resultWithDeletions.append(base)

    # Do insertions.
    resultWithDeletionsAndInsertions = list(resultWithDeletions)
    if DEBUG:
        debug(f'  resultWithDeletionsAndInsertions = '
              f'{"".join(resultWithDeletionsAndInsertions)}')
    for offset, bases in sorted(insertions.items()):
        if bases.count < minCoverage:
            base = lowCoverage
        else:
            base = leastAmbiguousFromCounts(bases.counts, threshold)

        deletionCount = sum(x <= offset for x in deletions)
        resultWithDeletionsAndInsertions.insert(offset - deletionCount, base)
        if DEBUG:
            debug(f'  resultWithDeletionsAndInsertions = '
                  f'{"".join(resultWithDeletionsAndInsertions)} inserted '
                  f'{base} before {offset - deletionCount}')

    if DEBUG:
        debug(f'orig                             = {"".join(orig)}')
        debug(f'result                           = {"".join(result)}')
        debug(f'resultWithDeletions              = '
              f'{"".join(resultWithDeletions)}')
        debug(f'resultWithDeletionsAndInsertions = '
              f'{"".join(resultWithDeletionsAndInsertions)}')
        debug(f'prefix                           = {"".join(prefix)}')
        debug(f'suffix                           = {"".join(suffix)}')

    return ''.join(prefix + resultWithDeletionsAndInsertions + suffix)


def addPairsInfo(pairs, query, qualities, referenceLength, correspondences,
                 deletions, insertions, insertionSet):
    """
    Add information about matched pairs of nucleotides.

    @param pairs: A C{list} of 2-C{tuple}s of query offset, reference offset.
        Either (but not both) member of each tuple might be C{None} to indicate
        an indel mismatch.
    @param query: A C{str} query DNA sequence.
    @param qualities: A C{list} of quality scores.
    @param correspondences: A C{defaultdict(list)}, to hold base, quality
        scores for when a query offset corresponds to a reference offset.
    @param deletions: A C{set} of C{int} reference offsets that are deleted in
        the query.
    @param insertions: A C{defaultdict(list)}, to hold base, quality
        scores for when a query contains an insertion to the reference.
    """
    # Find the offset of the sequence for the first member of the pair.
    leadingNoneCount = 0
    for _, referenceOffset in pairs:
        if referenceOffset is None:
            leadingNoneCount += 1
        else:
            break
    else:
        raise ValueError(
            'This is impossible in real life (though it can be triggered by a '
            'test), because no nucleotide in the read matches any location in '
            'the reference. I.e., the read did not map to the reference, '
            'which is impossible seeing as we are only working with reads '
            'that mapped.')

    if referenceOffset == 0:
        # The previous pairs correspond to nucleotides that match before
        # the start of the reference. So the actual offset (relative to
        # reference position 0) is negative.
        actualReferenceOffset = -leadingNoneCount
    else:
        actualReferenceOffset = referenceOffset

    insertCount = 0

    for queryOffset, referenceOffset in pairs:

        if queryOffset is None:
            # The query is missing something that is in the reference. So this
            # is a deletion from the reference.
            assert referenceOffset is not None
            # Sanity check.
            if referenceOffset is not None:
                assert referenceOffset == actualReferenceOffset
            deletions.add(actualReferenceOffset)
            actualReferenceOffset += 1

        elif (referenceOffset is None and
              actualReferenceOffset >= 0 and
              actualReferenceOffset < referenceLength):
            # The query has something that is not in the reference. So this
            # is an insertion to the reference.
            assert queryOffset is not None
            base = query[queryOffset]
            quality = qualities[queryOffset]
            offset = actualReferenceOffset + insertCount
            # insertCount = sum(insert <= offset for insert in insertionSet)
            insertCount = insertionSet.countLE(offset)
            offset = actualReferenceOffset + insertCount
            if DEBUG:
                debug(f'INSERT: {queryOffset=} {actualReferenceOffset=} '
                      f'{insertCount=}')
            insertions[offset].add(base, quality)
            insertionSet.add(offset)

        else:
            base = query[queryOffset]
            quality = qualities[queryOffset]
            correspondences[actualReferenceOffset].add(base, quality)
            actualReferenceOffset += 1

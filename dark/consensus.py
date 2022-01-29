import sys
from collections import defaultdict
from contextlib import contextmanager
import progressbar

from dark.dna import Bases
from dark.sam import (
    samfile, samReferences, UnequalReferenceLengthError,
    UnknownReference, UnspecifiedReference)

DEBUG = False


def debug(*msg):
    print(*msg, file=sys.stderr)


class Insertion:
    """
    Manage a collection of bases and quality scores for an insertion that
    starts at given genome site.

    @param insertionOffset: The C{int} offset of the insertion.
    """
    def __init__(self, insertionOffset):
        self.insertionOffset = insertionOffset
        self.anchorOffsets = []
        self.bases = []

    def __str__(self):
        result = []

        for anchorOffset, bases in zip(self.anchorOffsets, self.bases):
            result.append(f'  {anchorOffset}: {bases}')

        return (f'<Insertion at {self.insertionOffset}:\n' +
                '\n'.join(result) + '>')

    __repr__ = __str__

    def start(self, anchorOffset):
        """
        Start a new insertion.

        @param anchorOffset: The C{int} offset where the insertion should begin
            (in the reference). Or C{None} if the insertion should be placed
            immediately before {self.insertionOffset}.
        """
        self.anchorOffsets.append(anchorOffset)
        self.bases.append([])

    def append(self, base, quality):
        """
        Append a base and quality score to the current list.

        @param base: A C{str} nucleotide base.
        @param quality: An C{int} nucleotide quality score.
        @raise IndexError: If C{start} has not already been called.
        """
        self.bases[-1].append((base, quality))

    def resolve(self):
        """
        Resolve the insertion at this offset.

        @return: A C{list} of Bases() instances corresponding to the insertion.
        """
        # The length of the insertion will be at least the length of the
        # longest sequence of (base, quality) pairs we've been given for
        # this offset.
        nInsertions = max(len(bases) for bases in self.bases)

        # And we have to potentially add so that if some reads need to be
        # anchored at a certain starting offset.
        maxExtra = 0
        for anchorOffset, bases in zip(self.anchorOffsets, self.bases):
            if anchorOffset is not None:
                extra = (self.insertionOffset + nInsertions - anchorOffset -
                         len(bases) - 1)
                if DEBUG:
                    debug(f'{extra} = {self.insertionOffset} + '
                          f'{nInsertions} - {anchorOffset} - {len(bases)} - 1')
                if extra > maxExtra:
                    maxExtra = extra

        nInsertions += maxExtra

        # Make an array of Bases instances and add the (base, quality) pairs
        # to each one.
        insertion = []
        for _ in range(nInsertions):
            insertion.append(Bases())

        for anchorOffset, bases in zip(self.anchorOffsets, self.bases):
            startOffset = (
                nInsertions - len(bases) if anchorOffset is None else 0)
            for offset, (base, quality) in enumerate(bases):
                insertion[startOffset + offset].append(base, quality)

        if DEBUG:
            debug('insertion')
            for i in insertion:
                debug('\t', i)

        return insertion


@contextmanager
def maybeProgressBar(show, maxValue, prefix=None):
    """
    A context manager to maybe show a progress bar.

    @param show: If C{True}, yield a progress bar, else a Class with an
        C{update} method that does nothing.
    @param maxValue: The C{int} number of tasks to show progress for.
    @param prefix: A C{str} prefix, to appear at the start of the progress bar.
    """
    if show:
        with progressbar.ProgressBar(max_value=maxValue, prefix=prefix) as bar:
            yield bar
    else:
        class Bar:
            update = staticmethod(lambda _: None)
        yield Bar


def basesToConsensus(offsetBases, otherBases, originalOffsets, reference,
                     referenceLength, threshold, minCoverage, lowCoverage,
                     noCoverage):
    """
    Convert a C{dict} of C{Bases} into a consensus.

    @param offsetBases: A C{defaultdict} mapping C{int} offsets to
        C{Bases} instances.
    @param otherBases: A C{dict} mapping C{int} offsets to C{str} symbols
        for the consensus for offsets that we have no reads for.
    @param originalOffsets: A C{dict} mapping C{int} adjusted offsets to their
        C{int} original offsets.
    @param reference: A C{Read} instance giving the reference sequence, or
        C{None} (in which case C{lowCoverage} and C{noCoverage} may not be
        'reference'.
    @param referenceLength: The C{int} length of the reference (note that we
        may not have the reference sequence but we can still be passed its
        length, as found in the BAM file header).
    @param threshold: A C{float} threshold, as for C{consensusFromBAM}.
    @param minCoverage: An C{int} minimum number of reads that must cover a
        site for a consensus base to be called. If fewer reads cover a
        site, the C{lowCoverage} value is used.
    @param lowCoverage: A C{str} indicating what base to use when
        0 < N < minCoverage reads cover a site.
    @param noCoverage: A C{str} indicating what base to use when
        no reads cover the site.
    @return: a C{str} consensus sequence.
    """
    if DEBUG:
        debug(f'HERE: {offsetBases=}')
        debug(f'HERE: {otherBases=}')
    result = []
    allOffsets = set(offsetBases) | set(otherBases) | {0}
    minOffset = min(allOffsets)
    maxOffset = max(allOffsets)

    debug(f'basesToConsensus offset limits: {minOffset} - {maxOffset}')

    # prefixLen = -minOffset if minOffset < 0 else 0
    # suffixLen = (maxOffset - referenceLength + 1
    #              if maxOffset >= referenceLength else 0)
    #
    # if DEBUG:
    #     debug(f'  prefix len {prefixLen} suffix len {suffixLen}')

    noCoverageStr = (reference.sequence if noCoverage == 'reference' else
                     noCoverage * referenceLength)
    lowCoverageStr = (reference.sequence if lowCoverage == 'reference' else
                      lowCoverage * referenceLength)

    for offset in range(minOffset, maxOffset + 1):
        try:
            originalOffset = originalOffsets[offset]
            lowCoverageBase = lowCoverageStr[originalOffset]
            noCoverageBase = noCoverageStr[originalOffset]
        except (IndexError, KeyError):
            lowCoverageBase = noCoverageBase = '?'

        # Don't use try/except KeyError here as offsetBases is a defaultdict.
        if offset in sorted(offsetBases):
            bases = offsetBases[offset]

            if bases is None:
                if DEBUG:
                    debug(f'Offset {offset} was deleted. Skipping.')
            else:
                result.append(bases.consensus(threshold, minCoverage,
                                              lowCoverageBase, noCoverageBase))
        else:
            assert offset in otherBases, (
                f'Offset {offset} not found in offsetBases or otherBases.')
            if otherBases[offset] is None:
                if DEBUG:
                    debug(f'Adjusted {offset=} was deleted.')
            else:
                result.append(otherBases[offset])

    return ''.join(result)


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
    insertions = {}

    nReads = bam.count(contig=referenceId)
    with maybeProgressBar(progress, nReads, prefix='Reads: ') as bar:
        for readCount, read in enumerate(bam.fetch(contig=referenceId),
                                         start=1):
            assert not read.is_unmapped

            addPairsInfo(
                read.get_aligned_pairs(), read.query_sequence,
                ([1] * len(read.query_sequence) if ignoreQuality else
                 read.query_qualities),
                referenceLength, correspondences, deletions, insertions)

            if 1 or DEBUG:
                debug(f'read id  : {read.query_name}')
                debug('query    :', read.query_sequence)
                debug(f'cigar    : {read.cigarstring}')
                debug(f'match    : {read.reference_start}')
                debug(f'Pairs    : {read.get_aligned_pairs()}')
                debug(f'  {deletions=}')

                debug('correspondences:')
                for offset, item in sorted(correspondences.items()):
                    debug(f'  {offset}: {item}')

                debug('insertions:')
                for offset in sorted(insertions):
                    debug(f'  {insertions[offset]}')

            bar.update(readCount)

    conflicts = set(deletions) & set(insertions)
    if conflicts:
        if DEBUG:
            debug('CONFLICTING OFFSETS:',
                  ', '.join(sorted(map(str, conflicts))))

    # Do deletions.
    for offset in deletions:
        correspondences[offset] = None

    if DEBUG:
        debug(f'There are {len(insertions)} insertions.')

    # Do insertions.
    allOffsets = set(correspondences) | {0, referenceLength - 1}
    minOffset = min(allOffsets)
    maxOffset = max(allOffsets)

    debug(f'Offset limits: {minOffset} - {maxOffset}')

    debug(f'Insertion count {len(insertions)}')
    maxInsertionLength = -1
    for insertionOffset, insertion in insertions.items():
        for bases in insertion.bases:
            if len(bases) > maxInsertionLength:
                maxInsertionLength = len(bases)
        if maxInsertionLength > 350:
            debug(f'LONG (>350) INSERTION {insertion.insertionOffset}: '
                  f'{maxInsertionLength}')

    consensusBases = defaultdict(lambda: Bases())
    otherBases = {}
    originalOffsets = {}

    # debug(f'{correspondences=}')
    debug('correspondences:')
    for offset, item in sorted(correspondences.items()):
        debug(f'  {offset}: {item}')

    debug('insertions:')
    for offset in sorted(insertions):
        debug(f'  {insertions[offset]}')

    with maybeProgressBar(progress, maxOffset - minOffset + 1,
                          prefix='Sites: ') as bar:
        insertCount = deletionCount = 0
        for barCount, offset in enumerate(range(minOffset, maxOffset + 1),
                                          start=1):

            debug('OFFSET:', offset)

            if offset in insertions:
                insertion = insertions[offset]
                assert offset == insertion.insertionOffset
                insertionBases = insertion.resolve()
                insertionLength = len(insertionBases)
                if 1 or DEBUG:
                    debug('INSERTION:', insertion)
                    debug(f'  Resolved {insertionLength=}')
                    # debug(f'{insertionBases=}')

                for bases in insertionBases:
                    adjustedOffset = offset + insertCount - deletionCount
                    consensusBases[adjustedOffset] += bases
                    insertCount += 1
                    originalOffsets[adjustedOffset] = offset

            adjustedOffset = offset + insertCount - deletionCount
            originalOffsets[adjustedOffset] = offset

            if adjustedOffset > 30500:
                debug(f'HIGH!!! {adjustedOffset}')

            if offset in correspondences:
                if correspondences[offset] is None:
                    # This was deleted!
                    if DEBUG:
                        debug(f'  CONFLICTING OFFSET {offset}. Could not add '
                              f'{correspondences[offset]}')
                    if adjustedOffset in otherBases:
                        assert otherBases[adjustedOffset] is None
                    else:
                        otherBases[adjustedOffset] = None
                        deletionCount += 1
                else:
                    if DEBUG:
                        debug(f'Offset {offset} Adding correspondences bases '
                              f'{correspondences[offset]} to '
                              f'{consensusBases[adjustedOffset]} at offset '
                              f'{adjustedOffset}')
                    consensusBases[adjustedOffset] += correspondences[offset]
            else:
                debug(f'Offset {offset} not in correspondences')

            if adjustedOffset not in consensusBases:
                # We don't have any information (from the reads) for this
                # offset. If it's present in otherBases, it must be because
                # it was a deletion.
                if adjustedOffset in otherBases:
                    assert offset in deletions
                else:
                    otherBases[adjustedOffset] = (
                        reference.sequence[offset] if noCoverage == 'reference'
                        else noCoverage)

            bar.update(barCount)

    return basesToConsensus(consensusBases, otherBases, originalOffsets,
                            reference, referenceLength, threshold, minCoverage,
                            lowCoverage, noCoverage)


def addPairsInfo(pairs, query, qualities, referenceLength, correspondences,
                 deletions, insertions):
    """
    Add information about matched pairs of nucleotides.

    @param pairs: A C{list} of 2-C{tuple}s of query offset, reference offset.
        Either (but not both) member of each tuple might be C{None} to indicate
        an indel mismatch.
    @param query: A C{str} query DNA sequence.
    @param qualities: A C{list} of quality scores.
    @param correspondences: A C{defaultdict(list)}, to hold (base, quality)
        scores for when a query offset corresponds to a reference offset.
    @param deletions: A C{set} of C{int} reference offsets that are deleted in
        the query.
    @param insertions: A C{defaultdict(list)}, to hold (base, quality)
        scores for when a query contains an insertion to the reference.
    """
    if DEBUG:
        debug(f'Pairs: {pairs}')
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

    inInsertion = False

    for queryOffset, referenceOffset in pairs:
        # Sanity check.
        if referenceOffset is not None:
            assert referenceOffset == actualReferenceOffset

        if queryOffset is None:
            # The query is missing something that is in the reference. So this
            # is a deletion from the reference.
            assert referenceOffset is not None
            deletions.add(actualReferenceOffset)
            actualReferenceOffset += 1
            inInsertion = False

        elif (referenceOffset is None and
              actualReferenceOffset >= 0 and
              actualReferenceOffset < referenceLength):
            # The query has something that is not in the reference. So this
            # is an insertion to the reference.
            assert queryOffset is not None
            base = query[queryOffset]
            quality = qualities[queryOffset]

            if not inInsertion:
                inInsertion = True
                if actualReferenceOffset not in insertions:
                    insertions[actualReferenceOffset] = Insertion(
                        actualReferenceOffset)
                if leadingNoneCount:
                    # There were some leading None reference offsets, so
                    # these inserted bases will need to be positioned
                    # relative to the first non-None reference offset.
                    insertions[actualReferenceOffset].start(None)
                else:
                    insertions[actualReferenceOffset].start(
                        actualReferenceOffset)

            if DEBUG:
                debug(f'INSERT: {queryOffset=} {actualReferenceOffset=} '
                      f'{actualReferenceOffset=}')

            insertions[actualReferenceOffset].append(base, quality)

        else:
            base = query[queryOffset]
            quality = qualities[queryOffset]
            correspondences[actualReferenceOffset].append(base, quality)
            actualReferenceOffset += 1
            inInsertion = False

    if DEBUG:
        for insertionOffset, theseInsertions in insertions.items():
            debug(f'Pairs added. Offset {insertionOffset} has '
                  f'{theseInsertions}.')

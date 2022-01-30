import os
import sys
from collections import defaultdict
from contextlib import contextmanager
import progressbar
from pysam import CINS, CDEL, CSOFT_CLIP

from dark.cigar import (
    cigarTuplesToOperations, softClippedOffset, insertionOffset)
from dark.dna import Bases
from dark.sam import (
    samfile, samReferences, UnequalReferenceLengthError,
    UnknownReference, UnspecifiedReference, CONSUMES_REFERENCE)
from dark.utils import pct

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

    def readCount(self):
        """
        How many reads have an insertion at this offset?

        @return: An C{int} read count.
        """
        return len(self.bases)

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
def maybeProgressBar(show, maxValue, prefix):
    """
    A context manager to maybe show a progress bar.

    @param show: If C{True}, yield a progress bar, else a Class with an
        C{update} method that does nothing.
    @param maxValue: The C{int} number of tasks to show progress for.
    @param prefix: A C{str} prefix, to appear at the start of the progress bar.
    """
    if show and os.isatty(2):
        widgets = [
            progressbar.SimpleProgress(format='%(value_s)s/%(max_value_s)s'),
            progressbar.Percentage(format=' %(percentage)3d%%'),
            ' ',
            progressbar.Bar(marker='\x1b[33m#\x1b[39m'),
            ' ',
            progressbar.Timer(format='Elapsed: %(elapsed)s'),
            ' ',
            progressbar.ETA(format='ETA: %(eta)8s'),
        ]
        with progressbar.ProgressBar(max_value=maxValue, widgets=widgets,
                                     prefix=prefix) as bar:
            yield bar
    else:
        class Bar:
            update = staticmethod(lambda _: None)
        yield Bar


def basesToConsensus(offsetBases, otherBases, originalOffsets, reference,
                     referenceLength, threshold, minCoverage, lowCoverage,
                     noCoverage, deletionSymbol, progress):
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
    @parm deletionSymbol: The C{str} to insert in the consensus when a deleted
        site is detected.
    @param progress: If C{True}, display a progress bar on standard error.
    @return: a C{str} consensus sequence.
    """
    if DEBUG:
        debug(f'HERE: {offsetBases=}')
        debug(f'HERE: {otherBases=}')
    result = []
    allOffsets = set(offsetBases) | set(otherBases) | {0}
    minOffset = min(allOffsets)
    maxOffset = max(allOffsets)

    # debug(f'basesToConsensus offset limits: {minOffset} - {maxOffset}')

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

    with maybeProgressBar(progress, maxOffset - minOffset + 1,
                          'Consensus: ') as bar:
        for barCount, offset in enumerate(range(minOffset, maxOffset + 1),
                                          start=1):
            try:
                originalOffset = originalOffsets[offset]
                lowCoverageBase = lowCoverageStr[originalOffset]
                noCoverageBase = noCoverageStr[originalOffset]
            except (IndexError, KeyError):
                lowCoverageBase = noCoverageBase = '?'

            # Don't use try/except KeyError here as offsetBases is a
            # defaultdict.
            if offset in sorted(offsetBases):
                bases = offsetBases[offset]

                if bases is None:
                    result.append(deletionSymbol)
                else:
                    result.append(bases.consensus(
                        threshold, minCoverage, lowCoverageBase,
                        noCoverageBase))
            else:
                assert offset in otherBases, (
                    f'Offset {offset} not found in offsetBases or otherBases.')
                if otherBases[offset] is None:
                    result.append(deletionSymbol)
                else:
                    result.append(otherBases[offset])

            bar.update(barCount)

    return ''.join(result)


def consensusFromBAM(
        bamFilename, referenceId=None, reference=None, threshold=0.8,
        minCoverage=1, lowCoverage='reference', noCoverage='reference',
        deletionSymbol='-', deletionThreshold=0.5, ignoreQuality=False,
        insertionCountThreshold=5, strategy='majority',
        includeSoftClipped=False, progress=False):
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
    @parm deletionSymbol: The C{str} to insert in the consensus when a deleted
        site is detected.
    @param deletionThreshold: If some reads have a deletion at a site and some
        do not, call the site as a deletion if C{float} the fraction of reads
        with the deletion is at least this value.
    @param insertionCountThreshold: The C{int} number of reads that must have
        an insertion at an offset in order for the insertion to be called in
        the consensus.
    @param ignoreQuality: If C{True}, ignore quality scores.
    @param strategy: A C{str} consensus-making strategy (cuurently must be
        'majority').
    @param progress: If C{True}, display a progress bar on standard error.
    @param includeSoftClipped: Include information from read bases that were
        marked as soft-clipped by the algorithm that made the BAM file.
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
                minCoverage, lowCoverage, noCoverage, deletionSymbol,
                deletionThreshold, ignoreQuality, insertionCountThreshold,
                includeSoftClipped, progress)
        else:
            raise ValueError(f'Unknown consensus strategy {strategy!r}.')


def _fetchConsensus(bam, referenceId, reference, referenceLength, threshold,
                    minCoverage, lowCoverage, noCoverage, deletionSymbol,
                    deletionThreshold, ignoreQuality, insertionCountThreshold,
                    includeSoftClipped, progress):
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
    @parm deletionSymbol: The C{str} to insert in the consensus when a deleted
        site is detected.
    @param deletionThreshold: If some reads have a deletion at a site and some
        do not, call the site as a deletion if C{float} the fraction of reads
        with the deletion is at least this value.
    @param insertionCountThreshold: The C{int} number of reads that must have
        an insertion at an offset in order for the insertion to be called in
        the consensus.
    @param ignoreQuality: If C{True}, ignore quality scores.
        to.
    @param includeSoftClipped: Include information from read bases that were
        marked as soft-clipped by the algorithm that made the BAM file.
    @return: A C{str} consensus sequence.

    """
    correspondences = defaultdict(lambda: Bases())
    deletions = defaultdict(int)
    insertions = {}

    nReads = bam.count(contig=referenceId)
    with maybeProgressBar(progress, nReads, 'Reads    : ') as bar:
        for readCount, read in enumerate(bam.fetch(contig=referenceId),
                                         start=1):
            assert not read.is_unmapped

            addPairsInfo(
                read.get_aligned_pairs(),
                list(cigarTuplesToOperations(read.cigartuples)),
                read.query_sequence,
                ([1] * len(read.query_sequence) if ignoreQuality else
                 read.query_qualities), referenceLength, includeSoftClipped,
                correspondences, deletions, insertions)

            if DEBUG:
                debug(f'read id     : {read.query_name}')
                debug('query       :', read.query_sequence)
                debug('query len   :', len(read.query_sequence))
                debug(f'cigar       : {read.cigarstring}')
                debug(f'match offset: {read.reference_start}')
                debug(f'Pairs       : {read.get_aligned_pairs()}')
                debug(f'Deletions   : {deletions}')

                debug('correspondences:')
                for offset, item in sorted(correspondences.items()):
                    debug(f'  {offset}: {item}')

                debug('insertions:')
                for offset in sorted(insertions):
                    debug(f'  {insertions[offset]}')

            bar.update(readCount)

    conflicts = set(deletions) & set(insertions)
    if conflicts:
        debug('CONFLICTING OFFSETS:', ', '.join(sorted(map(str, conflicts))))

    _pileupCorrespondences = pileupCorrespondences(
        bam, referenceId, referenceLength, includeSoftClipped, progress)

    compareCorrespondences(correspondences, _pileupCorrespondences,
                           threshold, minCoverage)

    # Do deletions.
    if DEBUG and deletions:
        debug(f'There are {len(deletions)} deletions:')

    for offset, deletionCount in sorted(deletions.items()):
        count = correspondences[offset].count
        if count == 0 or deletionCount / count >= deletionThreshold:
            correspondences[offset] = None
            debug(f'  Offset {offset:5d}: deleted in {deletionCount:4d} '
                  f'present in {count:4d}')

    if DEBUG:
        debug('correspondences after deletions:')
        for offset, item in sorted(correspondences.items()):
            debug(f'  {offset}: {item}')

    # Do insertions.

    if DEBUG:
        debug(f'There are {len(insertions)} insertions.')

    allOffsets = set(correspondences) | {0, referenceLength - 1}
    minOffset = min(allOffsets)
    maxOffset = max(allOffsets)

    if DEBUG:
        debug(f'Offset limits: {minOffset} - {maxOffset}')
        debug(f'Insertion count {len(insertions)}')

    consensusBases = defaultdict(lambda: Bases())
    otherBases = {}
    originalOffsets = {}

    if 1 or DEBUG:
        debug('final insertions:')
        for offset in sorted(insertions):
            wanted = insertions[offset].readCount() >= insertionCountThreshold
            debug(f'  {"" if wanted else "un"}wanted {insertions[offset]}')

    with maybeProgressBar(progress, maxOffset - minOffset + 1,
                          'Collect  : ') as bar:
        insertCount = deletionCount = 0
        for barCount, offset in enumerate(range(minOffset, maxOffset + 1),
                                          start=1):

            # debug('OFFSET:', offset)

            if offset in insertions:
                insertion = insertions[offset]
                assert offset == insertion.insertionOffset
                if insertion.readCount() >= insertionCountThreshold:
                    insertionBases = insertion.resolve()
                    insertionLength = len(insertionBases)
                    if DEBUG:
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
                if DEBUG:
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
                            lowCoverage, noCoverage, deletionSymbol, progress)


def addPairsInfo(pairs, cigarOperations, query, qualities, referenceLength,
                 includeSoftClipped, correspondences, deletions, insertions):
    """
    Add information about matched pairs of nucleotides.

    @param pairs: A C{list} of 2-C{tuple}s of query offset, reference offset.
        Either (but not both) member of each tuple might be C{None} to indicate
        an indel mismatch.
    @param cigarOperations: A C{list} of CIGAR operations corresponding to the
        information in C{pairs}.
    @param query: A C{str} query DNA sequence.
    @param qualities: A C{list} of quality scores.
    @param includeSoftClipped: Include information from read bases that were
        marked as soft-clipped by the algorithm that made the BAM file.
    @param correspondences: A C{defaultdict(list)}, to hold (base, quality)
        scores for when a query offset corresponds to a reference offset.
    @param deletions: A C{set} of C{int} reference offsets that are deleted in
        the query.
    @param insertions: A C{defaultdict(list)}, to hold (base, quality)
        scores for when a query contains an insertion to the reference.
    """
    assert len(pairs) == len(cigarOperations)
    assert not any(pair == (None, None) for pair in pairs)

    # If the first reference offset is None, the first operation must be a
    # soft clip (not an insertion).
    # if pairs[0][1] is None:
    #     assert cigarOperations[0] == CSOFT_CLIP, (
    #         f'First CIGAR operation is {cigarOperations[0]}, '
    #         f'not a soft clip')

    inInsertion = False

    if DEBUG:
        debug(f'{pairs=}')

    for count, ((queryOffset, referenceOffset), cigarOperation) in enumerate(
            zip(pairs, cigarOperations)):

        if queryOffset is None:
            # The query is missing something that is in the reference. So this
            # is a deletion from the reference.
            assert cigarOperation == CDEL
            assert referenceOffset is not None
            deletions[referenceOffset] += 1
            inInsertion = False

        elif referenceOffset is None:
            base = query[queryOffset]
            quality = qualities[queryOffset]

            if cigarOperation == CINS:
                # The query has an insertion (relative to the reference).

                # A CIGAR string shouldn't start with an insertion, IMO.
                # Rather, in such a case, it must start with unmatched
                # (soft-clipped) bases.
                # assert lastReferenceOffset is not None

                lookedBack, iOffset = insertionOffset(count, pairs,
                                                      cigarOperations)
                if not inInsertion:
                    inInsertion = True
                    if iOffset not in insertions:
                        insertions[iOffset] = Insertion(iOffset)
                    insertions[iOffset].start(iOffset if lookedBack else None)

                insertions[iOffset].append(base, quality)
            else:
                assert cigarOperation == CSOFT_CLIP
                inInsertion = False
                if includeSoftClipped:
                    correspondences[
                        softClippedOffset(count, pairs, cigarOperations)
                    ].append(base, quality)
        else:
            # Query and reference offsets are both non-None.
            assert cigarOperation in CONSUMES_REFERENCE
            inInsertion = False
            base = query[queryOffset]
            quality = qualities[queryOffset]
            correspondences[referenceOffset].append(base, quality)


def pileupCorrespondences(bam, referenceId, referenceLength,
                          includeSoftClipped, progress):
    """
    Get bases for offsets that correspond to the reference using pysam's
    pileup funciton.

    @param bam: An open BAM file.
    @param referenceId: A C{str} reference name indicating which reference to
        find reads against in the BAM file.
    @param referenceLength: The C{int} length of the reference.
    @param includeSoftClipped: Include information from read bases that were
        marked as soft-clipped by the algorithm that made the BAM file.
    @param progress: If C{True}, display a progress bar on standard error.
    @return: A C{dict} keyed by C{int} reference offset, with C{Bases}
        instances as values.
    """
    correspondences = defaultdict(lambda: Bases())

    with maybeProgressBar(progress, referenceLength, 'Pileup   : ') as bar:
        for barCount, site in enumerate(bam.pileup(contig=referenceId),
                                        start=1):
            referenceOffset = site.reference_pos
            for read in site.pileups:
                if not includeSoftClipped:
                    cigarOperations = list(cigarTuplesToOperations(
                        read.alignment.cigartuples))
                queryOffset = read.query_position
                if queryOffset is None or (
                        not includeSoftClipped and
                        cigarOperations[queryOffset] == CSOFT_CLIP):
                    continue
                    # Make sure an empty Bases instance is placed in our
                    # results so there's a sign that the reference offset
                    # was returned by bam.pileup.
                    correspondences[referenceOffset]
                else:
                    alignment = read.alignment
                    correspondences[referenceOffset].append(
                        alignment.query_sequence[queryOffset],
                        alignment.query_qualities[queryOffset])

            bar.update(barCount)

    return correspondences


def compareCorrespondences(fetch, pileup, threshold, minCoverage):
    """
    Compare bases at corresponding offsets as obtained from the pysam fetch
    and pileup functions.

    @param fetch: A defaultdict(lambda: Bases()) keyed by C{int} offset.
    @param pileup: A defaultdict(lambda: Bases()) keyed by C{int} offset.
    @param threshold: A C{float} threshold, as in C{consensusFromBAM}.
    @param minCoverage: An C{int} minimum number of reads that must cover a
        site for a consensus base to be called.
    """
    fetchOffsets = set(offset for offset in fetch if offset >= 0)
    pileupOffsets = set(pileup)

    unmatched = pileupOffsets - fetchOffsets
    if unmatched:
        print(f'There are {len(unmatched)} sites obtained from pileup not '
              f'present in fetch:',
              file=sys.stderr)

        for offset in sorted(unmatched):
            print(f'  {offset + 1}: {pileup[offset]}', file=sys.stderr)

    unmatched = fetchOffsets - pileupOffsets
    if unmatched:
        print(f'There are {len(unmatched)} sites obtained from fetch not '
              f'present in pileup:', file=sys.stderr)

        for offset in sorted(unmatched):
            print(f'  {offset + 1}: {fetch[offset]}', file=sys.stderr)

    common = fetchOffsets & pileupOffsets
    if common:
        print(f'There are {len(common)} sites obtained from both fetch and '
              f'pileup.', file=sys.stderr)
        differenceCount = 0
        for offset in sorted(common):
            fetchBase = fetch[offset].consensus(
                threshold, minCoverage, 'LOW-COVERAGE', 'NO-COVERAGE')
            pileupBase = pileup[offset].consensus(
                threshold, minCoverage, 'LOW-COVERAGE', 'NO-COVERAGE')
            if fetchBase != pileupBase:
                differenceCount += 1
                print(f'  {offset + 1}: {fetchBase} != {pileupBase}',
                      file=sys.stderr)
                if fetch[offset]:
                    print(f'    Fetch : {fetch[offset]}', file=sys.stderr)
                if pileup[offset]:
                    print(f'    Pileup: {pileup[offset]}', file=sys.stderr)

        if differenceCount:
            print(f'There were differences called at '
                  f'{pct(differenceCount, len(common))} sites and '
                  f'{pct(len(common) - differenceCount, len(common))} were '
                  f'called identically.', file=sys.stderr)
        else:
            print(f'All {len(common)} sites were called identically.',
                  file=sys.stderr)

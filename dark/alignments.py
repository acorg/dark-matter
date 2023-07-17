from __future__ import annotations

import re
from typing import Union, Optional, Literal, Type, Generator, Callable
from dark.taxonomy import LineageFetcher
from dark.filter import TitleFilter
from dark.hsp import HSP, LSP
from dark.reads import Read, Reads
from dark.score import HigherIsBetterScore, LowerIsBetterScore


class Alignment:
    """
    Hold information about a read alignment.

    @param subjectLength: The C{int} length of the sequence a read matched
        against.
    @param subjectTitle: The C{str} title of the sequence a read matched
        against.
    """

    def __init__(self, subjectLength: int, subjectTitle: str):
        self.subjectLength = subjectLength
        self.subjectTitle = subjectTitle
        self.hsps: list[Union[HSP, LSP]] = []

    def addHsp(self, hsp: Union[HSP, LSP]) -> None:
        """
        Add an HSP (or LSP) to the list of HSPs for this alignment.

        @param hsp: A L{dark.hsp} (or subclass) instance.
        """
        self.hsps.append(hsp)


class ReadAlignments(list):
    """
    Holds information about the alignments for a read.

    @param read: A C{Read} instance.
    @param alignments: A C{list} of L{dark.alignment.Alignment} instances or
        C{None} if the read has no alignments.
    """

    def __init__(self, read: Read, alignments: Optional[list[Alignment]] = None):
        list.__init__(self)
        self.read = read
        if alignments:
            self.extend(alignments)


def bestAlignment(readAlignments: ReadAlignments) -> Alignment:
    """
    Find the best alignment for a read. This is the one whose first HSP
    has the best score.

    Note that the comparison of HSP score values is taken care of by
    the HSP class. This works whether higher or lower scores are
    considered better.

    @param readAlignments: a C{ReadAlignments} instance.
    @return: The alignment with the best first HSP score.

    """
    return max(readAlignments, key=lambda alignment: alignment.hsps[0])


class ReadsAlignmentsParams:
    """
    Holds information about how a ReadsAlignments instance was created.

    @param application: The C{str} name of the application that created
        the ReadsAlignments data.
    @param applicationParams: A C{dict} holding options and their values
        given to the application that created the ReadsAlignments data.
    @param subjectIsNucleotides: A C{bool} that indicates whether subject
        sequences (that are aligned against) are at the nucleotide level.
        If C{False}, the subject is assumed to be amino acids.
    @param scoreTitle: A C{str} to describe a score. E.g, "Bit value"
        or "E value".
    """

    def __init__(
        self,
        application: str,
        applicationParams: Optional[dict[str, str]] = None,
        subjectIsNucleotides: bool = True,
        scoreTitle: str = "Score",
    ):
        self.application = application
        self.applicationParams = applicationParams
        self.subjectIsNucleotides = subjectIsNucleotides
        self.scoreTitle = scoreTitle


class ReadsAlignmentsFilter:
    """
    Provide a filter for C{ReadsAlignments} instances.

    @param limit: An C{int} limit on the number of records to read.
    @param maxAlignmentsPerRead: An C{int} limit on the number of alignments a
        read may have in order not to be filtered out. Reads with a greater
        number of alignments will be elided. Pass 0 to filter out reads that
        did not match (i.e., align to) any subjects. Use C{None} for no
        max alignments filtering.
    @param minSequenceLen: Sequences of lesser length will be elided.
    @param maxSequenceLen: Sequences of greater length will be elided.
    @param minStart: HSPs that start before this offset in the matched
        sequence should not be returned.
    @param maxStop: HSPs that end after this offset in the matched sequence
        should not be returned.
    @param oneAlignmentPerRead: If C{True}, only keep the best
        alignment for each read.
    @param maxHspsPerHit: The maximum number of HSPs to keep for each
        alignment for each read.
    @param scoreCutoff: A C{float} score. Matches with scores that are not
        better than this score will be ignored.
    @param percentageIdenticalCutoff: A C{float} percentage identity (i.e.,
        0.0 to 100.0 NOT a fraction). Matches with percentage identity less
        than this will be ignored.
    @param percentagePositiveCutoff: A C{float} percentage (i.e.,
        0.0 to 100.0 NOT a fraction). Matches with a percentage positive less
        than this will be ignored. An AA match is considered positive (by
        DIAMOND) if its entry in the substitution (scoring) matrix is positive.
    @param whitelist: If not C{None}, a set of exact C{str} titles that are always
        acceptable (though the match info for a whitelist title may rule it
        out for other reasons).
    @param blacklist: If not C{None}, a set of exact C{str} titles that are never
        acceptable.
    @param whitelistFile: If not C{None}, a C{str} filename containing lines
        that give exact ids that are always acceptable.
    @param blacklistFile: If not C{None}, a C{str} filename containing lines
        that give exact ids that are never acceptable.
    @param titleRegex: A regex that sequence titles must match.
    @param negativeTitleRegex: A regex that sequence titles must not match.
    @param truncateTitlesAfter: A string that titles will be truncated
        beyond. If a truncated title has already been seen, that title will
        be elided.
    @param taxonomy: Either a C{str} name or an C{int} id of the taxonomic
        group on which should be filtered. eg 'Vira' will filter on
        viruses, while 11118 will filter on Coronaviridae.
    @param readIdRegex: A case-sensitive regex C{str} that read ids must
        match.
    @return: C{self}.
    """

    def __init__(
        self,
        limit: Optional[int] = None,
        maxAlignmentsPerRead: Optional[int] = None,
        minSequenceLen: Optional[int] = None,
        maxSequenceLen: Optional[int] = None,
        minStart: Optional[int] = None,
        maxStop: Optional[int] = None,
        oneAlignmentPerRead: bool = False,
        maxHspsPerHit: Optional[int] = None,
        scoreCutoff: Optional[float] = None,
        percentageIdenticalCutoff: Optional[float] = None,
        percentagePositiveCutoff: Optional[float] = None,
        whitelist: Optional[set[str]] = None,
        blacklist: Optional[set[str]] = None,
        whitelistFile: Optional[str] = None,
        blacklistFile: Optional[str] = None,
        titleRegex: Optional[str] = None,
        negativeTitleRegex: Optional[str] = None,
        truncateTitlesAfter: Optional[str] = None,
        taxonomy: Optional[Union[str, int]] = None,
        readIdRegex: Optional[str] = None,
    ):
        self.limit = limit
        self.maxAlignmentsPerRead = maxAlignmentsPerRead
        self.minSequenceLen = minSequenceLen
        self.maxSequenceLen = maxSequenceLen
        self.minStart = minStart
        self.maxStop = maxStop
        self.oneAlignmentPerRead = oneAlignmentPerRead
        self.maxHspsPerHit = maxHspsPerHit
        self.scoreCutoff = scoreCutoff
        self.percentageIdenticalCutoff = percentageIdenticalCutoff
        self.percentagePositiveCutoff = percentagePositiveCutoff
        self.titleFilter: Optional[TitleFilter]
        self.lineageFetcher: Optional[LineageFetcher]

        # If we've been asked to filter on matched sequence titles in any way,
        # build a title filter.
        if (
            whitelist
            or blacklist
            or whitelistFile
            or blacklistFile
            or titleRegex
            or negativeTitleRegex
            or truncateTitlesAfter
        ):
            self.titleFilter = TitleFilter(
                whitelist=whitelist,
                blacklist=blacklist,
                whitelistFile=whitelistFile,
                blacklistFile=blacklistFile,
                positiveRegex=titleRegex,
                negativeRegex=negativeTitleRegex,
                truncateAfter=truncateTitlesAfter,
            )
        else:
            self.titleFilter = None

        self.taxonomy = taxonomy
        self.lineageFetcher = None if taxonomy is None else LineageFetcher()
        self.readIdRegex = None if readIdRegex is None else re.compile(readIdRegex)
        self.count = 0

    def filter(self, readAlignments) -> Union[ReadAlignments, Literal[False]]:
        """
        Filter a read's alignments.

        @param readAlignments: A C{ReadAlignments} instance.
        @return: A C{ReadAlignments} instance if the passed
            C{readAlignments} is not filtered out, else C{False}.
        """

        # Implementation notes:
        #
        # 1. The order in which we carry out the filtering actions can make
        #    a big difference in the result of this function. The current
        #    ordering is based on what seems reasonable - it may not be the
        #    best way to do things. E.g., if maxHspsPerHit is 1 and there
        #    is a title regex, which should we perform first?
        #
        #    We perform filtering based on alignment before that based on
        #    HSPs. That's because there's no point filtering all HSPs for
        #    an alignment that we end up throwing away anyhow.
        #
        # 2. This function could be made faster if it first looked at its
        #    arguments and dynamically created an acceptance function
        #    (taking a readAlignments as an argument). The acceptance
        #    function would run without examining the desired filtering
        #    settings on each call the way the current code does.
        #
        # 3. A better approach with readIdRegex would be to allow the
        #    passing of a regex object. Then the caller would make the
        #    regex with whatever flags they liked (e.g., case insensitive).

        #
        # Alignment-only (i.e., non-HSP based) filtering.
        #
        if self.limit is not None and self.count == self.limit:
            return False

        # Does the read have too many alignments?
        if (
            self.maxAlignmentsPerRead is not None
            and len(readAlignments) > self.maxAlignmentsPerRead
        ):
            return False

        # Filter on the read id.
        if self.readIdRegex and self.readIdRegex.search(readAlignments.read.id) is None:
            return False

        if self.titleFilter:
            # Remove alignments against sequences whose titles are
            # unacceptable.
            wantedAlignments = []
            for alignment in readAlignments:
                if (
                    self.titleFilter.accept(alignment.subjectTitle)
                    != TitleFilter.REJECT
                ):
                    wantedAlignments.append(alignment)
            if wantedAlignments:
                readAlignments[:] = wantedAlignments
            else:
                return False

        # Only return alignments that are against sequences of the
        # desired length.
        minSequenceLen = self.minSequenceLen
        maxSequenceLen = self.maxSequenceLen
        if minSequenceLen is not None or maxSequenceLen is not None:
            wantedAlignments = []
            for alignment in readAlignments:
                length = alignment.subjectLength
                if not (
                    (minSequenceLen is not None and length < minSequenceLen)
                    or (maxSequenceLen is not None and length > self.maxSequenceLen)
                ):
                    wantedAlignments.append(alignment)
            if wantedAlignments:
                readAlignments[:] = wantedAlignments
            else:
                return False

        if self.taxonomy is not None:
            wantedAlignments = []
            for alignment in readAlignments:
                assert self.lineageFetcher
                lineage = self.lineageFetcher.lineage(alignment.subjectTitle)
                if lineage:
                    for taxonomyIdAndScientificName in lineage:
                        if self.taxonomy in taxonomyIdAndScientificName:
                            wantedAlignments.append(alignment)
                else:
                    # No lineage info was found. Keep the alignment
                    # since we can't rule it out.  We could add another
                    # option to control this.
                    wantedAlignments.append(alignment)
            if wantedAlignments:
                readAlignments[:] = wantedAlignments
            else:
                return False

        if self.oneAlignmentPerRead and readAlignments:
            readAlignments[:] = [bestAlignment(readAlignments)]

        #
        # From here on we do only HSP-based filtering.
        #

        # Throw out any unwanted HSPs due to maxHspsPerHit.
        if self.maxHspsPerHit is not None:
            for alignment in readAlignments:
                hsps = alignment.hsps
                if len(hsps) > self.maxHspsPerHit:
                    alignment.hsps = hsps[: self.maxHspsPerHit]

        # Throw out HSPs whose scores are not good enough.
        if self.scoreCutoff is not None:
            wantedAlignments = []
            for alignment in readAlignments:
                hsps = alignment.hsps
                wantedHsps = []
                for hsp in hsps:
                    if hsp.betterThan(self.scoreCutoff):
                        wantedHsps.append(hsp)
                if wantedHsps:
                    alignment.hsps = wantedHsps
                    wantedAlignments.append(alignment)
            if wantedAlignments:
                readAlignments[:] = wantedAlignments
            else:
                return False

        # Throw out HSPs whose percentage identical is not good enough.
        #
        # Note that if percentIdentical is None in the HSP it means that
        # we did not run DIAMOND or BLAST in a way that generated that
        # output, so we cannot eliminate such an HSP.
        if self.percentageIdenticalCutoff is not None:
            piCutoff = self.percentageIdenticalCutoff
            wantedAlignments = []
            for alignment in readAlignments:
                hsps = alignment.hsps
                wantedHsps = []
                for hsp in hsps:
                    pi = hsp.percentIdentical
                    if pi is None or pi >= piCutoff:
                        wantedHsps.append(hsp)
                if wantedHsps:
                    alignment.hsps = wantedHsps
                    wantedAlignments.append(alignment)
            if wantedAlignments:
                readAlignments[:] = wantedAlignments
            else:
                return False

        # Throw out HSPs whose percentage positive is not good enough.
        #
        # Note that if percentPositive is None in the HSP it means that
        # we did not run DIAMOND or BLAST in a way that generated that
        # output, so we cannot eliminate such an HSP.
        if self.percentagePositiveCutoff is not None:
            ppCutoff = self.percentagePositiveCutoff
            wantedAlignments = []
            for alignment in readAlignments:
                hsps = alignment.hsps
                wantedHsps = []
                for hsp in hsps:
                    pp = hsp.percentPositive
                    if pp is None or pp >= ppCutoff:
                        wantedHsps.append(hsp)
                if wantedHsps:
                    alignment.hsps = wantedHsps
                    wantedAlignments.append(alignment)
            if wantedAlignments:
                readAlignments[:] = wantedAlignments
            else:
                return False

        # Throw out HSPs that don't match in the desired place on the
        # matched sequence.
        minStart = self.minStart
        maxStop = self.maxStop
        if minStart is not None or maxStop is not None:
            wantedAlignments = []
            for alignment in readAlignments:
                hsps = alignment.hsps
                wantedHsps = []
                for hsp in hsps:
                    if not (
                        (minStart is not None and hsp.readStartInSubject < minStart)
                        or (maxStop is not None and hsp.readEndInSubject > maxStop)
                    ):
                        wantedHsps.append(hsp)
                if wantedHsps:
                    alignment.hsps = wantedHsps
                    wantedAlignments.append(alignment)
            if wantedAlignments:
                readAlignments[:] = wantedAlignments
            else:
                return False

        self.count += 1
        return readAlignments

    def close(self) -> None:
        """
        Close our lineage fetcher, if any.
        """
        if self.lineageFetcher:
            self.lineageFetcher.close()


class ReadsAlignments:
    """
    Provide for filtering for a collection of reads and their alignments.

    You probably will not use this class directly. Instead, write a
    subclass that implements iter.

    See L{dark.blast.alignments.BlastReadsAlignments} and
    L{dark.diamond.alignments.DiamondReadsAlignments} for examples.

    @param reads: A L{Reads} instance, containing the reads (query sequences)
        given to the application to create these matches.
    @param params: An instance of C{ReadsAlignmentsParams}, containing
        the details of the application that created the alignments.
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}, for comparing bit scores. If you
        are using e.g., BLAST e-values, pass LowerIsBetterScore instead.
    """

    def __init__(
        self,
        reads: Reads,
        params: ReadsAlignmentsParams,
        scoreClass: Union[
            Type[HigherIsBetterScore], Type[LowerIsBetterScore]
        ] = HigherIsBetterScore,
    ):
        self.reads = reads
        self.params = params
        self.scoreClass = scoreClass
        self._filters: list[Callable] = []

    def getSubjectSequence(self, title: str):
        """
        Obtain information about a sequence given its title.

        Must be implemented by a subclass, e.g., see
        L{blast.alignments.BlastReadsAlignments}.

        @param title: A C{str} sequence title from a BLAST or DIAMOND (etc.)
            match. Usually of the form
            'gi|63148399|gb|DQ011818.1| Description...'.
        @raise NotImplementedError: This method must be implemented by a
            subclass.
        """
        raise NotImplementedError(
            "getSubjectSequence must be implemented by " "a subclass"
        )

    def hsps(self) -> Generator[Union[HSP, LSP], None, None]:
        """
        Provide access to all HSPs for all alignments of all reads.

        @return: A generator that yields HSPs (or LSPs).
        """
        for readAlignments in self:
            for readAlignment in readAlignments:
                for hsp in readAlignment.hsps:
                    yield hsp

    def __iter__(self) -> Generator[ReadAlignments, None, None]:
        """
        Iterate through all readsAlignments, yielding the readAlignments that pass any
        filters that have been added.

        @return: A generator that yields readsAlignments.
        """
        for readAlignments in self.iter():
            for filterFunc in self._filters:
                filteredReadAlignments = filterFunc(readAlignments)
                if filteredReadAlignments is False:
                    break
                else:
                    readAlignments = filteredReadAlignments
            else:
                yield readAlignments

    def iter(self) -> Union[Generator[ReadAlignments, None, None], list]:
        """
        Placeholder to allow subclasses to provide readsAlignments.

        These might be extracted from a file. E.g., the
        C{dark.blast.alignments.blast.BlastReadsAlignments} class (a subclass
        of C{ReadsAlignments}) overrides this method to provide readsAlignments
        from files.

        @return: An iterable of C{ReadsAlignments} instances.
        """
        return []

    def filter(self, **kwargs) -> ReadsAlignments:
        """
        Add a filter to this C{readsAlignments}.

        @param kwargs: Keyword arguments, as accepted by
            C{ReadsAlignmentsFilter}.
        @return: C{self}
        """
        self._filters.append(ReadsAlignmentsFilter(**kwargs).filter)
        return self

from dark.taxonomy import LineageFetcher
from dark.filter import TitleFilter
from dark.score import HigherIsBetterScore


def bestAlignment(readAlignments):
    """
    Find the best alignment for a read. This is the one whose first HSP
    has the best score.

    Note that the comparison of HSP score values is taken care of by
    the HSP class. This works whether higher or lower scores are
    considered better.

    @param readAlignments: a C{ReadAlignments} instance.
    @return: The alignment with the best first HSP score.

    """
    return max(readAlignments.alignments,
               key=lambda alignment: alignment.hsps[0])


class Alignment(object):
    """
    Hold information about a read alignment.

    @param subjectLength: The C{int} length of the sequence the read hit
        against.
    @param subjectTitle: The C{str} title of the sequence the read hit against.
    """

    def __init__(self, subjectLength, subjectTitle):
        self.subjectLength = subjectLength
        self.subjectTitle = subjectTitle
        self.hsps = []

    def addHsp(self, hsp):
        """
        Add an HSP to the list of HSPs for this alignment.

        @param hsp: A L{dark.hsp} (or subclass) instance.
        """
        self.hsps.append(hsp)


class ReadAlignments(object):
    """
    Holds information about the alignments for a read.

    @param read: A C{Read} instance.
    @param alignments: A C{list} of L{dark.alignment.Alignment} instances.
    """
    def __init__(self, read, alignments):
        self.read = read
        self.alignments = alignments


class ReadsAlignmentsParams(object):
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

    def __init__(self, application, applicationParams=None,
                 subjectIsNucleotides=True, scoreTitle='Score'):
        self.application = application
        self.applicationParams = applicationParams
        self.subjectIsNucleotides = subjectIsNucleotides
        self.scoreTitle = scoreTitle


class ReadsAlignments(object):
    """
    Provide for filtering for a collection of reads and their alignments.

    You probably will not use this class directly. Instead, write a
    subclass that implements __iter__.
    See L{blast.alignments.BlastReadsAlignments} for an example.

    @param reads: A L{Reads} instance, containing the reads (sequences)
        given to the application to create these hits.
    @param params: An instance of C{ReadsAlignmentsParams}, containing
        the details of the application that created the alignments.
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}, for comparing bit scores. If you
        are using e.g., BLAST e-values, pass LowerIsBetterScore instead.
    """

    def __init__(self, reads, params, scoreClass=HigherIsBetterScore):
        self.reads = reads
        self.params = params
        self.scoreClass = scoreClass

    def __iter__(self):
        """
        Must be implemented by a subclass, e.g., see
        L{blast.alignments.BlastReadsAlignments}.
        """
        raise NotImplementedError('__iter__ must be implemented by a subclass')

    def getSequence(self, title):
        """
        Obtain information about a sequence given its title.

        Must be implemented by a subclass, e.g., see
        L{blast.alignments.BlastReadsAlignments}.

        @param title: A C{str} sequence title from a BLAST hit. Of the form
            'gi|63148399|gb|DQ011818.1| Description...'.
        @return: A C{SeqIO.read} instance.
        """
        raise NotImplementedError('getSequence must be implemented by a '
                                  'subclass')

    def hsps(self):
        """
        Provide access to all HSPs for all alignments of all reads.

        @return: A generator that yields HSPs (or LSPs).
        """
        for readAlignment in self:
            for alignment in readAlignment.alignments:
                for hsp in alignment.hsps:
                    yield hsp

    def filter(self, limit=None, minSequenceLen=None, maxSequenceLen=None,
               minStart=None, maxStop=None,
               oneAlignmentPerRead=False, maxHspsPerHit=None,
               scoreCutoff=None, whitelist=None, blacklist=None,
               titleRegex=None, negativeTitleRegex=None,
               truncateTitlesAfter=None, taxonomy=None):
        """
        @param limit: An C{int} limit on the number of records to read.
        @param minSequenceLen: sequences of lesser length will be elided.
        @param maxSequenceLen: sequences of greater length will be elided.
        @param minStart: HSPs that start before this offset in the hit sequence
            should not be returned.
        @param maxStop: HSPs that end after this offset in the hit sequence
            should not be returned.
        @param oneAlignmentPerRead: if C{True}, only keep the best
            alignment for each read.
        @param maxHspsPerHit: The maximum number of HSPs to keep for each
            alignment for each read.
        @param scoreCutoff: A C{float} score. Hits with scores that are not
            better than this score will be ignored.
        @param whitelist: If not C{None}, a set of exact titles that are always
            acceptable (though the hit info for a whitelist title may rule it
            out for other reasons).
        @param blacklist: If not C{None}, a set of exact titles that are never
            acceptable.
        @param titleRegex: a regex that sequence titles must match.
        @param negativeTitleRegex: a regex that sequence titles must not match.
        @param truncateTitlesAfter: specify a string that titles will be
            truncated beyond. If a truncated title has already been seen, that
            title will be elided.
        @param taxonomy: a C{str} of the taxonomic group on which should be
            filtered. eg 'Vira' will filter on viruses.
        @return: A generator that yields L{Hit} instances.
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
        #    (taking a hit as an argument). The acceptance function would
        #    run without examining the above arguments for each hit the way
        #    the current code does.

        #
        # Alignment-only (i.e., non-HSP based) filtering.
        #

        # If we've been asked to filter on hit sequence titles in any way,
        # build a title filter.
        if (whitelist or blacklist or titleRegex or negativeTitleRegex or
                truncateTitlesAfter):
            titleFilter = TitleFilter(
                whitelist=whitelist, blacklist=blacklist,
                positiveRegex=titleRegex, negativeRegex=negativeTitleRegex,
                truncateAfter=truncateTitlesAfter)
        else:
            titleFilter = None

        if taxonomy:
            lineageFetcher = LineageFetcher()

        count = 0
        for hit in self:
            if limit is not None and count == limit:
                return

            if titleFilter:
                # Remove alignments against sequences whose titles are
                # unacceptable.
                wantedAlignments = []
                for alignment in hit.alignments:
                    if (titleFilter.accept(alignment.subjectTitle) !=
                            TitleFilter.REJECT):
                        wantedAlignments.append(alignment)
                if wantedAlignments:
                    hit.alignments = wantedAlignments
                else:
                    continue

            # Only return alignments that are against sequences of the
            # desired length.
            if minSequenceLen is not None or maxSequenceLen is not None:
                wantedAlignments = []
                for alignment in hit.alignments:
                    length = alignment.subjectLength
                    if not ((minSequenceLen is not None and
                             length < minSequenceLen) or
                            (maxSequenceLen is not None and
                             length > maxSequenceLen)):
                        wantedAlignments.append(alignment)
                if wantedAlignments:
                    hit.alignments = wantedAlignments
                else:
                    continue

            if taxonomy:
                wantedAlignments = []
                for alignment in hit.alignments:
                    lineage = lineageFetcher.lineage(alignment.subjectTitle)
                    if lineage:
                        if taxonomy in lineage:
                            wantedAlignments.append(alignment)
                    else:
                        # No lineage info was found. Keep the alignment
                        # since we can't rule it out.  We could add another
                        # option to control this.
                        wantedAlignments.append(alignment)
                if wantedAlignments:
                    hit.alignments = wantedAlignments
                else:
                    continue

            if oneAlignmentPerRead and hit.alignments:
                hit.alignments = [bestAlignment(hit)]

            #
            # From here on we do only HSP-based filtering.
            #

            # Throw out any unwanted HSPs due to maxHspsPerHit.
            if maxHspsPerHit is not None:
                for alignment in hit.alignments:
                    hsps = alignment.hsps
                    if len(hsps) > maxHspsPerHit:
                        alignment.hsps = hsps[:maxHspsPerHit]

            # Throw out HSPs whose scores are not good enough.
            if scoreCutoff is not None:
                wantedAlignments = []
                for alignment in hit.alignments:
                    hsps = alignment.hsps
                    wantedHsps = []
                    for hsp in hsps:
                        if hsp.betterThan(scoreCutoff):
                            wantedHsps.append(hsp)
                    if wantedHsps:
                        alignment.hsps = wantedHsps
                        wantedAlignments.append(alignment)
                if wantedAlignments:
                    hit.alignments = wantedAlignments
                else:
                    continue

            # Throw out HSPs that don't match in the desired place on the
            # hit sequence.
            if minStart is not None or maxStop is not None:
                wantedAlignments = []
                for alignment in hit.alignments:
                    hsps = alignment.hsps
                    wantedHsps = []
                    for hsp in hsps:
                        if not ((minStart is not None and
                                 hsp.readStartInSubject < minStart)
                                or (maxStop is not None and
                                    hsp.readEndInSubject > maxStop)):
                            wantedHsps.append(hsp)
                    if wantedHsps:
                        alignment.hsps = wantedHsps
                        wantedAlignments.append(alignment)
                if wantedAlignments:
                    hit.alignments = wantedAlignments
                else:
                    continue

            yield hit
            count += 1

        if taxonomy:
            lineageFetcher.close()

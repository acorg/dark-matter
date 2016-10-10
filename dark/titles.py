from collections import defaultdict, Counter

from dark.reads import Reads
from dark.utils import median
from dark.filter import ReadSetFilter
from dark.intervals import ReadIntervals


def titleCounts(readsAlignments):
    """
    Count the number of times each title in a readsAlignments instance is
    matched. This is useful for rapidly discovering what titles were matched
    and with what frequency.

    @param readsAlignments: A L{dark.alignments.ReadsAlignments} instance.
    @return: A C{dict} whose keys are titles and whose values are the integer
        counts of the number of reads that matched that title.
    """

    titles = defaultdict(int)
    for readAlignments in readsAlignments:
        for alignment in readAlignments:
            titles[alignment.subjectTitle] += 1
    return titles


class TitleAlignment(object):
    """
    Hold information about a read's HSPs for a title alignment.

    @param read: The C{Read} that aligned.
    @param hsps: A C{list} of L{dark.hsp.HSP} (or subclass) instance.
    """

    def __init__(self, read, hsps):
        self.read = read
        self.hsps = hsps


class TitleAlignments(list):
    """
    Holds information about a list of alignments against a sequence.

    @param subjectTitle: The C{str} title of the sequence the read matched
        against.
    @param subjectLength: The C{int} length of the sequence the read matched
        against.
    """

    def __init__(self, subjectTitle, subjectLength):
        # TODO: Do we need the title in here?
        self.subjectTitle = subjectTitle
        self.subjectLength = subjectLength

    def addAlignment(self, alignment):
        """
        Add an alignment to the list of alignments that matched this title.

        @param alignment: A L{TitleAlignment} instance.
        """
        self.append(alignment)

    def reads(self):
        """
        Find the set of reads matching a this title.

        @return: An instance of C{dark.reads.Reads}.
        """
        reads = Reads()
        for alignment in self:
            reads.add(alignment.read)
        return reads

    def readCount(self):
        """
        Find out how many reads aligned to this title.

        @return: The C{int} number of reads that aligned to this title.
        """
        return len(self)

    def hspCount(self):
        """
        How many HSPs were there in total for all the alignments to a title.

        @return: The C{int} number of HSPs for the alignments to this title.
        """
        return sum(len(alignment.hsps) for alignment in self)

    def readIds(self):
        """
        Find the set of read ids that matched the title.

        @return: A C{set} of read ids that aligned to this title.
        """
        return set(alignment.read.id for alignment in self)

    def hsps(self):
        """
        Get all HSPs for the alignments to a title.

        @return: A generator yielding L{dark.hsp.HSP} instances.
        """
        return (hsp for alignment in self for hsp in alignment.hsps)

    def bestHsp(self):
        """
        Find the HSP with the best score.

        @raise ValueError: If there are no HSPs.
        @return: The C{dark.hsp.HSP} instance (or a subclass) with the best
        score.
        """
        return max(hsp for hsp in self.hsps())

    def worstHsp(self):
        """
        Find the HSP with the worst score.

        @raise ValueError: If there are no HSPs.
        @return: The C{dark.hsp.HSP} instance (or a subclass) with the worst
        score.
        """
        return min(hsp for hsp in self.hsps())

    def hasScoreBetterThan(self, score):
        """
        Is there an HSP with a score better than a given value?

        @return: A C{bool}, C{True} if there is at least one HSP in the
        alignments for this title with a score better than C{score}.
        """
        # Note: Do not assume that HSPs in an alignment are sorted in
        # decreasing order (as they are in BLAST output). If we could
        # assume that, we could just check the first HSP in each alignment.
        for hsp in self.hsps():
            if hsp.betterThan(score):
                return True
        return False

    def medianScore(self):
        """
        Find the median score for the HSPs in the alignments that match
        this title.

        @raise ValueError: If there are no HSPs.
        @return: The C{float} median score of HSPs in alignments matching the
            title.
        """
        return median([hsp.score.score for hsp in self.hsps()])

    def coverage(self):
        """
        Get the fraction of this title sequence that is matched by its reads.

        @return: The C{float} fraction of the title sequence matched by its
            reads.
        """
        intervals = ReadIntervals(self.subjectLength)
        for hsp in self.hsps():
            intervals.add(hsp.subjectStart, hsp.subjectEnd)
        return intervals.coverage()

    def residueCounts(self, convertCaseTo='upper'):
        """
        Count residue frequencies at all sequence locations matched by reads.

        @param convertCaseTo: A C{str}, 'upper', 'lower', or 'none'.
            If 'none', case will not be converted (both the upper and lower
            case string of a residue will be present in the result if they are
            present in the read - usually due to low complexity masking).
        @return: A C{dict} whose keys are C{int} offsets into the title
            sequence and whose values are C{Counters} with the residue as keys
            and the count of that residue at that location as values.
        """
        if convertCaseTo == 'none':
            convert = lambda x: x
        elif convertCaseTo == 'lower':
            convert = str.lower
        elif convertCaseTo == 'upper':
            convert = str.upper
        else:
            raise ValueError(
                "convertCaseTo must be one of 'none', 'lower', or 'upper'")

        counts = defaultdict(Counter)

        for titleAlignment in self:
            read = titleAlignment.read
            for hsp in titleAlignment.hsps:
                for (subjectOffset, residue, inMatch) in read.walkHSP(hsp):
                    counts[subjectOffset][convert(residue)] += 1

        return counts


class TitlesAlignments(dict):
    """
    Holds (as a dictionary) a set of titles, each with its alignments.

    @param readsAlignments: A L{dark.alignments.ReadsAlignments} instance.
    @param scoreClass: A class to hold and compare scores. If C{None},
        the score class from readsAlignments will be used.
    @param readSetFilter: An instance of dark.filter.ReadSetFilter, or C{None}.
        This can be used to pass a previously used title filter for ongoing
        use in filtering.
    @param importReadsAlignmentsTitles: If C{True}, titles from
        C{readsAlignments} will be added to self. This argument is only used
        by the filtering function to make a new instance without reading its
        titles.
    """

    def __init__(self, readsAlignments, scoreClass=None, readSetFilter=None,
                 importReadsAlignmentsTitles=True):
        dict.__init__(self)
        self.readsAlignments = readsAlignments
        self.scoreClass = scoreClass or readsAlignments.scoreClass
        self.readSetFilter = readSetFilter
        if importReadsAlignmentsTitles:
            for readAlignments in readsAlignments:
                for alignment in readAlignments:
                    title = alignment.subjectTitle
                    try:
                        titleAlignments = self[title]
                    except KeyError:
                        titleAlignments = self[title] = TitleAlignments(
                            title, alignment.subjectLength)
                    titleAlignments.addAlignment(
                        TitleAlignment(readAlignments.read, alignment.hsps))

    def addTitle(self, title, titleAlignments):
        """
        Add a new title to self.

        @param title: A C{str} title.
        @param titleAlignments: An instance of L{TitleAlignments}.
        @raises KeyError: If the title is already present.
        """
        if title in self:
            raise KeyError('Title %r already present in TitlesAlignments '
                           'instance.' % title)
        else:
            self[title] = titleAlignments

    def filter(self, minMatchingReads=None, minMedianScore=None,
               withScoreBetterThan=None, minNewReads=None, minCoverage=None,
               maxTitles=None, sortOn='maxScore'):
        """
        Filter the titles in self to create another TitlesAlignments.

        @param minMatchingReads: titles that are matched by fewer reads
            are unacceptable.
        @param minMedianScore: sequences that are matched with a median
            bit score that is less are unacceptable.
        @param withScoreBetterThan: if the best score for a title is not
            as good as this value, the title is not acceptable.
        @param minNewReads: The C{float} fraction of its reads by which a new
            title's read set must differ from the read sets of all previously
            seen titles in order for this title to be considered acceptably
            different (and therefore interesting).
        @param minCoverage: The C{float} minimum fraction of the title sequence
            that must be matched by at least one read.
        @param maxTitles: A non-negative C{int} maximum number of titles to
            keep. If more titles than this are present, titles will be sorted
            (according to C{sortOn}) and only the best will be retained.
        @param sortOn: A C{str} attribute to sort on, used only if C{maxTitles}
            is not C{None}. See the C{sortTitles} method below for the legal
            values.
        @raise: C{ValueError} if C{maxTitles} is less than zero or the value of
            C{sortOn} is unknown.
        @return: A new L{TitlesAlignments} instance containing only the
            matching titles.
        """
        # Use a ReadSetFilter only if we're checking that read sets are
        # sufficiently new.
        if minNewReads is None:
            readSetFilter = None
        else:
            if self.readSetFilter is None:
                self.readSetFilter = ReadSetFilter(minNewReads)
            readSetFilter = self.readSetFilter

        result = TitlesAlignments(
            self.readsAlignments, self.scoreClass, self.readSetFilter,
            importReadsAlignmentsTitles=False)

        if maxTitles is not None and len(self) > maxTitles:
            if maxTitles < 0:
                raise ValueError('maxTitles (%r) cannot be negative.' %
                                 maxTitles)
            else:
                # There are too many titles. Make a sorted list of them so
                # we loop through them (below) in the desired order and can
                # break when/if we've reached the maximum. We can't just
                # take the first maxTitles titles from the sorted list now,
                # as some of those titles might later be discarded by the
                # filter and then we'd return a result with fewer titles
                # than we should.
                titles = self.sortTitles(sortOn)
        else:
            titles = self.keys()

        for title in titles:
            # Test max titles up front, as it may be zero.
            if maxTitles is not None and len(result) == maxTitles:
                break

            titleAlignments = self[title]
            if (minMatchingReads is not None and
                    titleAlignments.readCount() < minMatchingReads):
                continue

            # To compare the median score with another score, we must
            # convert both values to instances of the score class used in
            # this data set so they can be compared without us needing to
            # know if numerically greater scores are considered better or
            # not.
            if (minMedianScore is not None and
                    self.scoreClass(titleAlignments.medianScore()) <
                    self.scoreClass(minMedianScore)):
                continue

            if (withScoreBetterThan is not None and not
                    titleAlignments.hasScoreBetterThan(withScoreBetterThan)):
                continue

            if (minCoverage is not None and
                    titleAlignments.coverage() < minCoverage):
                continue

            if (readSetFilter and not
                    readSetFilter.accept(title, titleAlignments)):
                continue

            result.addTitle(title, titleAlignments)

        return result

    def hsps(self):
        """
        Get all HSPs for all the alignments for all titles.

        @return: A generator yielding L{dark.hsp.HSP} instances.
        """
        return (hsp for titleAlignments in self.values()
                for alignment in titleAlignments for hsp in alignment.hsps)

    def sortTitles(self, by):
        """
        Sort titles by a given attribute and then by title.

        @param by: A C{str}, one of 'length', 'maxScore', 'medianScore',
            'readCount', or 'title'.
        @return: A sorted C{list} of titles.
        """
        # First sort titles by the secondary key, which is always the title.
        titles = sorted(iter(self))

        # Then sort on the primary key (if any).
        if by == 'length':
            return sorted(
                titles,  reverse=True,
                key=lambda title: self[title].subjectLength)
        if by == 'maxScore':
            return sorted(
                titles, reverse=True, key=lambda title: self[title].bestHsp())
        if by == 'medianScore':
            return sorted(
                titles, reverse=True,
                key=lambda title: self.scoreClass(self[title].medianScore()))
        if by == 'readCount':
            return sorted(
                titles, reverse=True,
                key=lambda title: self[title].readCount())
        if by == 'title':
            return titles

        raise ValueError('Sort attribute must be one of "length", "maxScore", '
                         '"medianScore", "readCount", "title".')

import numpy as np
from collections import defaultdict

from dark.reads import Reads
from dark.filter import ReadSetFilter


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

        @return: The C{dark.hsp.HSP} instance (or a subclass) with the best
        score.
        """
        return max(hsp for hsp in self.hsps())

    def worstHsp(self):
        """
        Find the HSP with the worst score.

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
        Find out the median score for the HSPs in the alignments that match
        this title.

        @return: The C{float} median score of HSPs in alignments matching the
            title.
        """
        return np.median([hsp.score.score for hsp in self.hsps()])


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
        by the filtering function to make an new instance without reading its
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
        @param titleAlignments: An instance of L{TitleAlignment}.
        @raises KeyError: If the title is already present.
        """
        if title in self:
            raise KeyError('Title %r already present in TitlesAlignments '
                           'instance.' % title)
        else:
            self[title] = titleAlignments

    def filter(self, minMatchingReads=None, minMedianScore=None,
               withScoreBetterThan=None, minNewReads=None):
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

        for title, titleAlignments in self.iteritems():
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
        return (hsp for titleAlignments in self.itervalues()
                for alignment in titleAlignments for hsp in alignment.hsps)

    def sortTitles(self, by):
        """
        Sort titles by a given attribute and then by title.

        @param by: A C{str}, one of 'length', 'maxScore', 'medianScore',
            'readCount', or 'title'.
        @return: A sorted C{list} of titles.
        """

        def makeScoreCmp(attr):
            """
            Create a sorting comparison function that sorts first in reverse
            on the result (as a score) of calling the passed attribute
            function and then in ascending order on title.
            """
            def compare(title1, title2):
                value1 = self.scoreClass(getattr(self[title1], attr)())
                value2 = self.scoreClass(getattr(self[title2], attr)())
                result = cmp(value2, value1)
                if result == 0:
                    result = cmp(title1, title2)
                return result
            return compare

        def makeHspCmp(attr):
            """
            Create a sorting comparison function that sorts first in reverse
            on the result (as a score) of calling the passed attribute
            function and then in ascending order on title.
            """
            def compare(title1, title2):
                value1 = getattr(self[title1], attr)()
                value2 = getattr(self[title2], attr)()
                result = cmp(value2, value1)
                if result == 0:
                    result = cmp(title1, title2)
                return result
            return compare

        def makeNumericCmp(attr):
            """
            Create a sorting comparison function that sorts first in reverse
            on the passed attribute (which may be a function or numeric) and
            then in ascending order on title.
            """
            def compare(title1, title2):
                value1 = getattr(self[title1], attr)
                value2 = getattr(self[title2], attr)
                if callable(value1) and callable(value2):
                    result = cmp(value2(), value1())
                else:
                    result = cmp(value2, value1)
                if result == 0:
                    result = cmp(title1, title2)
                return result
            return compare

        if by == 'length':
            return sorted(self.iterkeys(), cmp=makeNumericCmp('subjectLength'))
        if by == 'maxScore':
            return sorted(self.iterkeys(), cmp=makeHspCmp('bestHsp'))
        if by == 'medianScore':
            return sorted(self.iterkeys(), cmp=makeScoreCmp('medianScore'))
        if by == 'readCount':
            return sorted(self.iterkeys(), cmp=makeNumericCmp('readCount'))
        if by == 'title':
            return sorted(self.iterkeys())

        raise ValueError('Sort attribute must be one of "length", "maxScore", '
                         '"medianScore", "readCount", "title".')

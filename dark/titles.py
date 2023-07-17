from collections import defaultdict, Counter

from dark.utils import median
from dark.filter import ReadSetFilter, TitleFilter
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


class TitleAlignment:
    """
    Hold information about a read's HSPs for a title alignment.

    @param read: The C{Read} that aligned.
    @param hsps: A C{list} of L{dark.hsp.HSP} (or subclass) instances.
    """

    def __init__(self, read, hsps):
        self.read = read
        self.hsps = hsps

    def toDict(self):
        """
        Get information about a title alignment as a dictionary.

        @return: A C{dict} representation of the title aligment.
        """
        return {
            "hsps": [hsp.toDict() for hsp in self.hsps],
            "read": self.read.toDict(),
        }


class TitleAlignments(list):
    """
    Holds information about a list of alignments against a sequence.

    @param subjectTitle: The C{str} title of the sequence the read matched
        against.
    @param subjectLength: The C{int} length of the sequence the read matched
        against.
    """

    def __init__(self, subjectTitle, subjectLength):
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
        Find the set of reads matching this title.

        @return: A generator that yields C{dark.reads.Read} instances (or one
            of its subclasses).
        """
        return (alignment.read for alignment in self)

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

    def coverageCounts(self):
        """
        For each location in the title sequence, return a count of how many
        times that location is covered by a read.
        """
        intervals = ReadIntervals(self.subjectLength)
        for hsp in self.hsps():
            intervals.add(hsp.subjectStart, hsp.subjectEnd)
        return intervals.coverageCounts()

    def coverageInfo(self):
        """
        Return information about the bases found at each location in our title
        sequence.

        @return: A C{dict} whose keys are C{int} subject offsets and whose
            values are unsorted lists of (score, base) 2-tuples, giving all the
            bases from reads that matched the subject at subject location,
            along with the bit score of the matching read.
        """
        result = defaultdict(list)

        for titleAlignment in self:
            for hsp in titleAlignment.hsps:
                score = hsp.score.score
                for subjectOffset, base, _ in titleAlignment.read.walkHSP(
                    hsp, includeWhiskers=False
                ):
                    result[subjectOffset].append((score, base))

        return result

    def residueCounts(self, convertCaseTo="upper"):
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
        if convertCaseTo == "none":

            def convert(x):
                return x

        elif convertCaseTo == "lower":
            convert = str.lower
        elif convertCaseTo == "upper":
            convert = str.upper
        else:
            raise ValueError("convertCaseTo must be one of 'none', 'lower', or 'upper'")

        counts = defaultdict(Counter)

        for titleAlignment in self:
            read = titleAlignment.read
            for hsp in titleAlignment.hsps:
                for subjectOffset, residue, inMatch in read.walkHSP(hsp):
                    counts[subjectOffset][convert(residue)] += 1

        return counts

    def summary(self):
        """
        Summarize the alignments for this subject.

        @return: A C{dict} with C{str} keys:
            bestScore: The C{float} best score of the matching reads.
            coverage: The C{float} fraction of the subject genome that is
                matched by at least one read.
            hspCount: The C{int} number of hsps that match the subject.
            medianScore: The C{float} median score of the matching reads.
            readCount: The C{int} number of reads that match the subject.
            subjectLength: The C{int} length of the subject.
            subjectTitle: The C{str} title of the subject.
        """
        return {
            "bestScore": self.bestHsp().score.score,
            "coverage": self.coverage(),
            "hspCount": self.hspCount(),
            "medianScore": self.medianScore(),
            "readCount": self.readCount(),
            "subjectLength": self.subjectLength,
            "subjectTitle": self.subjectTitle,
        }

    def toDict(self):
        """
        Get information about the title's alignments as a dictionary.

        @return: A C{dict} representation of the title's aligments.
        """
        return {
            "titleAlignments": [titleAlignment.toDict() for titleAlignment in self],
            "subjectTitle": self.subjectTitle,
            "subjectLength": self.subjectLength,
        }


class TitlesAlignments(dict):
    """
    Holds (as a dictionary) a set of titles, each with its alignments.

    @param readsAlignments: A L{dark.alignments.ReadsAlignments} instance.
    @param scoreClass: A class to hold and compare scores. If C{None},
        the score class from readsAlignments will be used.
    @param readSetFilter: An instance of dark.filter.ReadSetFilter, or C{None}.
        This can be used to pass a previously used title filter for ongoing
        use in filtering.
    """

    def __init__(self, readsAlignments, scoreClass=None, readSetFilter=None):
        dict.__init__(self)
        self.readsAlignments = readsAlignments
        self.scoreClass = scoreClass or readsAlignments.scoreClass
        self.readSetFilter = readSetFilter
        for readAlignments in readsAlignments:
            for alignment in readAlignments:
                title = alignment.subjectTitle
                try:
                    titleAlignments = self[title]
                except KeyError:
                    titleAlignments = self[title] = TitleAlignments(
                        title, alignment.subjectLength
                    )
                titleAlignments.addAlignment(
                    TitleAlignment(readAlignments.read, alignment.hsps)
                )

    def addTitle(self, title, titleAlignments):
        """
        Add a new title to self.

        @param title: A C{str} title.
        @param titleAlignments: An instance of L{TitleAlignments}.
        @raises KeyError: If the title is already present.
        """
        if title in self:
            raise KeyError(
                "Title %r already present in TitlesAlignments " "instance." % title
            )
        else:
            self[title] = titleAlignments

    def filter(
        self,
        minMatchingReads=None,
        maxMatchingReads=None,
        minMedianScore=None,
        withScoreBetterThan=None,
        minNewReads=None,
        minCoverage=None,
        maxTitles=None,
        sortOn="maxScore",
        titleRegex=None,
        negativeTitleRegex=None,
    ):
        """
        Filter the titles in self.

        @param minMatchingReads: titles that are matched by fewer reads
            are unacceptable.
        @param maxMatchingReads: titles that are matched by more reads
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
        @param titleRegex: A regex that read ids must match.
        @param negativeTitleRegex: A regex that read ids must not match.
        @raise: C{ValueError} if C{maxTitles} is less than zero or the value of
            C{sortOn} is unknown.
        @return: C{self}, with non-matching titles removed.
        """
        # Use a ReadSetFilter only if we're checking that read sets are
        # sufficiently new.
        if minNewReads is None:
            readSetFilter = None
        else:
            if self.readSetFilter is None:
                self.readSetFilter = ReadSetFilter(minNewReads)
            readSetFilter = self.readSetFilter

        if maxTitles is not None and len(self) > maxTitles:
            if maxTitles < 0:
                raise ValueError("maxTitles (%r) cannot be negative." % maxTitles)
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
            titles = list(self)

        if titleRegex or negativeTitleRegex:
            titleFilter = TitleFilter(
                positiveRegex=titleRegex, negativeRegex=negativeTitleRegex
            )
        else:
            titleFilter = None

        titlesToKeep = set()

        for title in titles:
            # Test max titles up front, as it may be zero.
            if maxTitles is not None and len(titlesToKeep) == maxTitles:
                break

            # Test positive and negative regexps.
            if titleFilter and titleFilter.accept(title) == TitleFilter.REJECT:
                continue

            titleAlignments = self[title]
            if (
                minMatchingReads is not None
                and titleAlignments.readCount() < minMatchingReads
            ):
                continue

            if (
                maxMatchingReads is not None
                and titleAlignments.readCount() > maxMatchingReads
            ):
                continue

            # To compare the median score with another score, we must
            # convert both values to instances of the score class used in
            # this data set so they can be compared without us needing to
            # know if numerically greater scores are considered better or
            # not.
            if minMedianScore is not None and self.scoreClass(
                titleAlignments.medianScore()
            ) < self.scoreClass(minMedianScore):
                continue

            if (
                withScoreBetterThan is not None
                and not titleAlignments.hasScoreBetterThan(withScoreBetterThan)
            ):
                continue

            if minCoverage is not None and titleAlignments.coverage() < minCoverage:
                continue

            if readSetFilter and not readSetFilter.accept(title, titleAlignments):
                continue

            titlesToKeep.add(title)

        # Filter self.
        for title in list(self):
            if title not in titlesToKeep:
                del self[title]

        return self

    def hsps(self):
        """
        Get all HSPs for all the alignments for all titles.

        @return: A generator yielding L{dark.hsp.HSP} instances.
        """
        return (
            hsp
            for titleAlignments in self.values()
            for alignment in titleAlignments
            for hsp in alignment.hsps
        )

    def sortTitles(self, by):
        """
        Sort titles by a given attribute and then by title.

        @param by: A C{str}, one of 'length', 'maxScore', 'medianScore',
            'readCount', or 'title'.
        @raise ValueError: If an unknown C{by} value is given.
        @return: A sorted C{list} of titles.
        """
        # First sort titles by the secondary key, which is always the title.
        titles = sorted(iter(self))

        # Then sort on the primary key (if any).
        if by == "length":
            return sorted(
                titles, reverse=True, key=lambda title: self[title].subjectLength
            )
        if by == "maxScore":
            return sorted(titles, reverse=True, key=lambda title: self[title].bestHsp())
        if by == "medianScore":
            return sorted(
                titles,
                reverse=True,
                key=lambda title: self.scoreClass(self[title].medianScore()),
            )
        if by == "readCount":
            return sorted(
                titles, reverse=True, key=lambda title: self[title].readCount()
            )
        if by == "title":
            return titles

        raise ValueError(
            'Sort attribute must be one of "length", "maxScore", '
            '"medianScore", "readCount", "title".'
        )

    def summary(self, sortOn=None):
        """
        Summarize all the alignments for this title.

        @param sortOn: A C{str} attribute to sort titles on. One of 'length',
            'maxScore', 'medianScore', 'readCount', or 'title'.
        @raise ValueError: If an unknown C{sortOn} value is given.
        @return: A generator that yields C{dict} instances as produced by
            C{TitleAlignments} (see class earlier in this file), sorted by
            C{sortOn}.
        """
        titles = self if sortOn is None else self.sortTitles(sortOn)

        for title in titles:
            yield self[title].summary()

    def tabSeparatedSummary(self, sortOn=None):
        """
        Summarize all the alignments for this title as multi-line string with
        TAB-separated values on each line.

        @param sortOn: A C{str} attribute to sort titles on. One of 'length',
            'maxScore', 'medianScore', 'readCount', or 'title'.
        @raise ValueError: If an unknown C{sortOn} value is given.
        @return: A newline-separated C{str}, each line with a summary of a
            title. Each summary line is TAB-separated.
        """
        # The order of the fields returned here is somewhat arbitrary. The
        # subject titles are last because they are so variable in length.
        # Putting them last makes it more likely that the initial columns in
        # printed output will be easier to read down.
        #
        # Note that post-processing scripts will be relying on the field
        # ordering here.  So you can't just add fields. It's probably safe
        # to add them at the end, but be careful / think.
        #
        # A TAB-separated file can easily be read by awk using e.g.,
        # awk 'BEGIN {FS = "\t"} ...'

        result = []
        for titleSummary in self.summary(sortOn):
            result.append(
                "\t".join(
                    [
                        "%(coverage)f",
                        "%(medianScore)f",
                        "%(bestScore)f",
                        "%(readCount)d",
                        "%(hspCount)d",
                        "%(subjectLength)d",
                        "%(subjectTitle)s",
                    ]
                )
                % titleSummary
            )
        return "\n".join(result)

    def toDict(self):
        """
        Get information about the titles alignments as a dictionary.

        @return: A C{dict} representation of the titles aligments.
        """
        return {
            "scoreClass": self.scoreClass.__name__,
            "titles": dict(
                (title, titleAlignments.toDict())
                for title, titleAlignments in self.items()
            ),
        }

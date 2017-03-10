from math import log
from collections import Counter


class ReadIntervals(object):
    """
    Hold information about the set of reads that match a subject.

    @param targetLength: The C{int} length of the target sequence that the
        reads are against.
    """

    EMPTY = 0
    FULL = 1

    def __init__(self, targetLength):
        self._targetLength = targetLength
        self._intervals = []

    def add(self, start, end):
        """
        Add the start and end offsets of a matching read.

        @param start: The C{int} start offset of the read match in the subject.
        @param end: The C{int} end offset of the read match in the subject.
            This is Python-style: the end offset is not included in the match.
        """
        assert start <= end
        self._intervals.append((start, end))

    def walk(self):
        """
        Get the non-overlapping read intervals that match the subject.

        @return: A generator that produces (TYPE, (START, END)) tuples, where
            where TYPE is either self.EMPTY or self.FULL and (START, STOP) is
            the interval. The endpoint (STOP) of the interval is not considered
            to be in the interval. I.e., the interval is really [START, STOP).
        """
        intervals = sorted(self._intervals)

        def nextFull():
            start, stop = intervals.pop(0)
            while intervals:
                if intervals[0][0] <= stop:
                    _, thisStop = intervals.pop(0)
                    if thisStop > stop:
                        stop = thisStop
                else:
                    break
            return (start, stop)

        if intervals:
            # If the first interval (read) starts after zero, yield an
            # initial empty section to get us to the first interval.
            if intervals[0][0] > 0:
                yield (self.EMPTY, (0, intervals[0][0]))

            while intervals:
                # Yield a full interval followed by an empty one (if there
                # is another interval pending).
                lastFull = nextFull()
                yield (self.FULL, lastFull)
                if intervals:
                    yield (self.EMPTY, (lastFull[1], intervals[0][0]))

            # Yield the final empty section, if any.
            if lastFull[1] < self._targetLength:
                yield (self.EMPTY, (lastFull[1], self._targetLength))

        else:
            yield (self.EMPTY, (0, self._targetLength))

    def coverage(self):
        """
        Get the fraction of a subject is matched by its set of reads.

        @return: The C{float} fraction of a subject matched by its reads.
        """
        if self._targetLength == 0:
            return 0.0

        coverage = 0
        for (intervalType, (start, end)) in self.walk():
            if intervalType == self.FULL:
                # Adjust start and end to ignore areas where the read falls
                # outside the target.
                coverage += (min(end, self._targetLength) - max(0, start))
        return float(coverage) / self._targetLength

    def coverageCounts(self):
        """
        For each location in the subject, return a count of how many times that
        location is covered by a read.

        @return: a C{Counter} where the keys are the C{int} locations on the
            subject and the value is the number of times the location is
            covered by a read.
        """
        coverageCounts = Counter()
        for start, end in self._intervals:
            coverageCounts.update(range(max(0, start),
                                        min(self._targetLength, end)))
        return coverageCounts


class OffsetAdjuster(object):
    """
    A class that knows how to adjust the offsets in a normalized HSP according
    to the overall set of reads being plotted.

    @param intervals: An instance of L{ReadIntervals}.
    @param base: The C{float} logarithmic base to use when adjusting empty
        spaces in the hit sequence.
    """

    def __init__(self, intervals=None, base=2.0):
        self._adjustments = []  # Pairs of (X offset, adjustment).
        if intervals:
            divisor = log(base)
            for (intervalType, (start, stop)) in intervals.walk():
                if intervalType == ReadIntervals.EMPTY:
                    width = stop - start
                    logWidth = log(width) / divisor
                    self._adjustments.append((stop, width - logWidth))

    def adjustments(self):
        """
        Provide the adjustment values for this instance.

        @return: A C{list} of (X offset, adjustment) values, where the X
            offset is an C{int} and the adjustment is a C{float}.
        """
        return self._adjustments

    def _reductionForOffset(self, offset):
        """
        Calculate the total reduction for a given X axis offset.

        @param offset: The C{int} offset.
        @return: The total C{float} reduction that should be made for this
            offset.
        """
        reduction = 0
        for (thisOffset, thisReduction) in self._adjustments:
            if offset >= thisOffset:
                reduction += thisReduction
            else:
                break
        return reduction

    def adjustOffset(self, offset):
        """
        Adjust a single X offset.

        @param offset: The C{int} offset to adjust.
        @return: The C{float} adjusted offset.
        """
        return offset - self._reductionForOffset(offset)

    def adjustHSP(self, hsp):
        """
        Adjust the read and subject start and end offsets in an HSP.

        @param hsp: a L{dark.hsp.HSP} or L{dark.hsp.LSP} instance.
        """
        reduction = self._reductionForOffset(
            min(hsp.readStartInSubject, hsp.subjectStart))

        hsp.readEndInSubject = hsp.readEndInSubject - reduction
        hsp.readStartInSubject = hsp.readStartInSubject - reduction
        hsp.subjectEnd = hsp.subjectEnd - reduction
        hsp.subjectStart = hsp.subjectStart - reduction

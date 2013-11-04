from math import log


class ReadIntervals(object):
    """

    targetLength: the length of the target sequence that these reads
    were against.
    """

    EMPTY = 0
    FULL = 1

    def __init__(self, targetLength):
        self._targetLength = targetLength
        self._intervals = []

    def add(self, start, end):
        assert start <= end
        self._intervals.append((start, end))

    def walk(self):
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


class OffsetAdjuster(object):
    """
    A class that knows how to adjust the offsets in a normalized HSP according
    to the overall set of reads being plotted.

    intervals: an instance of ReadIntervals.
    base: the logarithmic base to use when adjusting empty spaces in the hit
        sequence.
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
        return self._adjustments

    def _reductionForOffset(self, offset):
        """Calculate the total reduction for a given X axis offset."""
        reduction = 0
        for (thisOffset, thisReduction) in self._adjustments:
            if offset >= thisOffset:
                reduction += thisReduction
            else:
                break
        return reduction

    def adjustOffset(self, offset):
        """Adjust a single X offset."""
        return offset - self._reductionForOffset(offset)

    def adjustNormalizedHSP(self, hsp):
        """
        Adjust the query and subject start and end offsets in a normalized HSP,
        returning a new normalized HSP.

        hsp: a normalized HSP, as produced by normalizeHSP in hsps.py
        """
        reduction = self._reductionForOffset(hsp['queryStart'])

        return {
            'queryEnd': hsp['queryEnd'] - reduction,
            'queryStart': hsp['queryStart'] - reduction,
            'subjectEnd': hsp['subjectEnd'] - reduction,
            'subjectStart': hsp['subjectStart'] - reduction,
        }

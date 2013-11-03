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

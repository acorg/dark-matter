from operator import itemgetter


class ColorsForCounts:
    """
    Maintain a collection of count thresholds and colors with methods to get a
    color or a CSS name for a count.

    @param colors: An C{iterable} of space separated "value color" strings,
        such as ["100 red", "200 rgb(23, 190, 207)", "700 #CF3CF3"]. Or C{None}
        if no colors (other than C{defaultColor}) should be used.
    @param defaultColor: The C{str} color to use for counts that do not reach
        the lowest count threshold for any color in C{colors}.
    @raise ValueError: If an incorrect count/color pair is found in C{colors}.
    """

    def __init__(self, colors, defaultColor="black"):
        thresholds = set()
        result = []
        if colors:
            for colorInfo in colors:
                fields = colorInfo.split(None, 1)
                if len(fields) == 2:
                    threshold, color = fields
                    try:
                        threshold = int(threshold)
                    except ValueError:
                        raise ValueError(
                            "color arguments must be given as space-separated "
                            'pairs of "count color" where the count is an '
                            "integer threshold. Your value (%r) was not "
                            "an integer." % threshold
                        )
                    if threshold < 0:
                        raise ValueError(
                            "color arguments must be given as space-separated "
                            'pairs of "count color" where the count is '
                            "non-negative. Your value (%r) is less than 0." % threshold
                        )
                    if threshold in thresholds:
                        raise ValueError(
                            "repeated color argument count (%d)." % threshold
                        )

                    result.append((threshold, color))
                    thresholds.add(threshold)
                else:
                    raise ValueError(
                        "color arguments must be given as space-separated "
                        'pairs of "value color". Your value (%r) does not '
                        "contain a space." % colorInfo
                    )

            result.sort(key=itemgetter(0), reverse=True)

        if not result or result[-1][0] > 0:
            result.append((0, defaultColor))

        self.colors = tuple(result)

    def thresholdToCssName(self, threshold):
        """
        Turn a count threshold into a string that can be used as a CSS
        class name.

        @param threshold: The C{int} threshold.
        @raise ValueError: If the threshold is not an C{int}.
        @return: A C{str} CSS class name.
        """
        return "threshold-%d" % threshold

    def thresholdForCount(self, count):
        """
        Get the best threshold for a specific count.

        @param count: An C{int} count.
        @return: The first C{int} threshold that the given count is at least
            as big as.
        """
        assert count >= 0, "Count (%d) cannot be negative." % count
        for threshold, _ in self.colors:
            if count >= threshold:
                return threshold

        raise ValueError("This should never happen! Last threshold is not 0?")

    def colorForCount(self, count):
        """
        Get the color for a count.

        @param count: An C{int} count.
        @return: The C{str} color for the count.
        """
        assert count >= 0, "Count (%d) cannot be negative." % count
        for threshold, color in self.colors:
            if count >= threshold:
                return color

        raise ValueError("This should never happen! Last threshold is not 0?")

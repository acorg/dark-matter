import re
from math import ceil
from collections import OrderedDict

from dark.simplify import simplifyTitle


class TitleFilter(object):
    """
    Provide an acceptance test for sequence titles.

    @param whitelist: If not C{None}, a set of exact titles that are always
        acceptable.
    @param blacklist: If not C{None}, a set of exact titles that are never
        acceptable.
    @param positiveRegex: If not C{None}, a regex that sequence titles must
        match.
    @param negativeRegex: If not C{None}, a regex that sequence titles must
        not match.
    @param truncateAfter: A C{str} that titles will be truncated beyond. If
        a truncated title has already been seen, that title will no longer
        be acceptable.
    """

    REJECT = 0
    WHITELIST_ACCEPT = 1
    DEFAULT_ACCEPT = 2

    def __init__(self, whitelist=None, blacklist=None, positiveRegex=None,
                 negativeRegex=None, truncateAfter=None):
        self._whitelist = whitelist
        self._blacklist = blacklist
        if truncateAfter is None:
            self._truncated = None
        else:
            self._truncateAfter = truncateAfter
            self._truncated = {}

        if positiveRegex is None:
            self._positiveRegex = None
        else:
            self._positiveRegex = re.compile(positiveRegex, re.I)

        if negativeRegex is None:
            self._negativeRegex = None
        else:
            self._negativeRegex = re.compile(negativeRegex, re.I)

    def accept(self, title):
        """
        Return a value (see below) to indicate if a title is acceptable (and,
        if so, in what way).

        @param title: A C{str} sequence title.
        @return: An C{int} to indicate an acceptable title or not. This will be

            C{self.REJECT} if the title is unacceptable.
            C{self.WHITELIST_ACCEPT} if the title is whitelisted.
            C{self.DEFAULT_ACCEPT} if the title is acceptable by default.

            These three values are needed so our caller can distinguish between
            the two reasons for acceptance.
        """
        if self._whitelist and title in self._whitelist:
            return self.WHITELIST_ACCEPT

        if self._blacklist and title in self._blacklist:
            return self.REJECT

        if self._positiveRegex and self._positiveRegex.search(title) is None:
            return self.REJECT

        if (self._negativeRegex and
                self._negativeRegex.search(title) is not None):
            return self.REJECT

        if self._truncated is not None:
            truncated = simplifyTitle(title, self._truncateAfter)
            if truncated in self._truncated:
                # We've already seen this (truncated) title. Reject unless
                # this is the original title that we truncated to make this
                # entry. That title must continue to be accepted.
                if self._truncated[truncated] == title:
                    return self.DEFAULT_ACCEPT
                else:
                    return self.REJECT
            else:
                self._truncated[truncated] = title

        return self.DEFAULT_ACCEPT


class ReadSetFilter(object):
    """
    Provide an acceptance test based on sequence read set.

    @param minNew: The C{float} fraction of its reads by which a new read set
        must differ from all previously seen read sets in order to be
        considered acceptably different.
    """

    def __init__(self, minNew):
        self._minNew = minNew
        # Use an OrderedDict so that each time we walk through self._titles
        # we do it in the same order. This makes our runs deterministic /
        # reproducible.
        self._titles = OrderedDict()

    def accept(self, title, titleAlignments):
        """
        Return C{True} if the read id set in C{titleAlignments} is sufficiently
        different from all previously seen read sets.

        @param title: A C{str} sequence title.
        @param titleAlignments: An instance of L{TitleAlignment}.
        @return: A C{bool} indicating whether a title has an acceptably novel
            read set or not.
        """

        # Sanity check: titles can only be passed once.
        assert title not in self._titles, (
            'Title %r seen multiple times.' % title)

        readIds = titleAlignments.readIds()
        newReadsRequired = ceil(self._minNew * len(readIds))

        for readSet, invalidatedTitles in self._titles.values():
            if len(readIds - readSet) < newReadsRequired:
                # Add this title to the set of titles invalidated by this
                # previously seen read set.
                invalidatedTitles.append(title)
                return False

        # Remember the new read set and an empty list of invalidated titles.
        self._titles[title] = (readIds, [])

        return True

    def invalidates(self, title):
        """
        Report on which other titles were invalidated by a given title.

        @param title: A C{str} sequence title.
        @return: A C{list} of titles that the passed title invalidated.
        """
        try:
            return self._titles[title][1]
        except KeyError:
            return []

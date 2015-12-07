from __future__ import division

import string
from os.path import basename


def numericallySortFilenames(names):
    """
    Sort (ascending) a list of file names by their numerical prefixes.
    The number sorted on is the numeric prefix of the basename of
    the given filename. E.g., '../output/1.json.bz2' will sort before
    '../output/10.json.bz2'.

    @param: A C{list} of file names, each of whose basename starts with a
        string of digits.
    @return: The sorted C{list} of full file names.
    """

    def numericPrefix(name):
        """
        Find any numeric prefix of C{name} and return it as an C{int}.

        @param: A C{str} file name, whose name possibly starts with digits.
        @return: The C{int} number at the start of the name, else 0 if
            there are no leading digits in the name.
        """
        count = 0
        for ch in name:
            if ch in string.digits:
                count += 1
            else:
                break
        return 0 if count == 0 else int(name[0:count])

    return sorted(names, key=lambda name: numericPrefix(basename(name)))


# Set up a median function that works on different Python versions and on pypy
# (even if numpypy is still broken).
#
# First try numpy, and check that an ndarray has a partition method (to avoid a
# shortcoming in pypy numpy, see
# https://bitbucket.org/pypy/numpy/issues/12/median-calling-ndarraypartition
# If numpy seems ok, use it. Else try for the Python 3.4 median function, else
# write our own.

import numpy as np
_a = np.array([1])

try:
    _a.partition
except AttributeError:
    # Cannot use numpy. Try for the new (Python 3.4) built-in function.
    try:
        from statistics import median as _median
    except ImportError:
        def _median(l):
            """
            Find the median of a set of numbers.

            @param l: A C{list} of numeric values.
            @return: The median value in C{l}.
            """
            l = sorted(l)
            n = len(l)
            if n % 2:
                return l[n // 2]
            else:
                n //= 2
                return (l[n - 1] + l[n]) / 2
else:
    # We have a working numpy. Just make sure we raise on an empty list as
    # numpy returns numpy.nan
    _median = np.median

del _a


def median(l):
    """
    Find the median of a set of numbers.

    @param l: A C{list} of numeric values.
    @raise ValueError: If C{l} is empty.
    @return: The median value in C{l}.
    """

    # Handle the zero length case ahead of time because various median
    # implementations do different things. numpy returns nan, whereas Python
    # 3.4 and later raise StatisticsError

    if len(l) == 0:
        # Match the empty-argument exception message of Python's built-in max
        # and min functions.
        raise ValueError('arg is an empty sequence')
    else:
        return _median(l)

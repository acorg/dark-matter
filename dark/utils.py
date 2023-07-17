import os
import re
import string
import bz2
import gzip
from pathlib import Path
from os.path import basename
from contextlib import contextmanager
from re import compile
from statistics import median as _median
from typing import List, Optional


def numericallySortFilenames(names: List[str]) -> List[str]:
    """
    Sort (ascending) a list of file names by their numerical prefixes.
    The number sorted on is the numeric prefix of the basename of
    the given filename. E.g., '../output/1.json.bz2' will sort before
    '../output/10.json.bz2'.

    @param: A C{list} of file names, each of whose basename starts with a
        string of digits.
    @return: The sorted C{list} of full file names.
    """

    def numericPrefix(name: str) -> int:
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


def median(numbers):
    """
    Find the median of a set of numbers.

    @param numbers: A C{list} of numeric values.
    @raise ValueError: If C{l} is empty.
    @return: The median of C{numbers}.
    """
    if len(numbers) == 0:
        # Match the empty-argument exception message of Python's built-in max
        # and min functions.
        raise ValueError("arg is an empty sequence")
    else:
        return _median(numbers)


@contextmanager
def asHandle(fileNameOrHandle, mode="rt", encoding="UTF-8"):
    """
    Decorator for file opening that makes it easy to open compressed files
    and which can be passed an already-open file handle or a file name.
    Based on L{Bio.File.as_handle}.

    @param fileNameOrHandle: Either a C{str} or a file handle.
    @param mode: The C{str} mode to use for opening the file.
    @param encoding: The C{str} encoding to use when opening the file.
    @return: A generator that can be turned into a context manager via
        L{contextlib.contextmanager}.
    """
    if isinstance(fileNameOrHandle, (Path, str)):
        fileNameOrHandle = str(fileNameOrHandle)
        if fileNameOrHandle.endswith(".gz"):
            yield gzip.open(fileNameOrHandle, mode=mode, encoding=encoding)
        elif fileNameOrHandle.endswith(".bz2"):
            yield bz2.open(fileNameOrHandle, mode=mode, encoding=encoding)
        else:
            # Putting mode=mode, encoding=encoding into the following
            # causes a hard-to-understand error from the mocking library.
            with open(fileNameOrHandle) as fp:
                yield fp
    else:
        yield fileNameOrHandle


_rangeRegex = compile(r"^\s*(\d+)(?:\s*-\s*(\d+))?\s*$")


# Note: the parseRangeExpression, which uses the following function as a
#       helper (see below) is more general / powerful than this function on
#       its own.
def parseRangeString(s, convertToZeroBased=False):
    """
    Parse a range string of the form 1-5,12,100-200.

    @param s: A C{str} specifiying a set of numbers, given in the form of
        comma separated numeric ranges or individual indices.
    @param convertToZeroBased: If C{True} all indices will have one
        subtracted from them.
    @raise ValueError: If the range in C{s} cannot be parsed.
    @return: A C{set} of all C{int}s in the specified set.
    """

    result = set()
    for _range in s.split(","):
        match = _rangeRegex.match(_range)
        if match:
            start, end = match.groups()
            start = int(start)
            if end is None:
                end = start
            else:
                end = int(end)
            if start > end:
                start, end = end, start
            if convertToZeroBased:
                result.update(range(start - 1, end))
            else:
                result.update(range(start, end + 1))
        else:
            raise ValueError(
                "Illegal range %r. Ranges must single numbers or "
                "number-number." % _range
            )

    return result


# The following matches range expressions such as
#
#   3
#   3-4
#   3-4,5
#   3-4,6,9-10,99
#
# including any embedded whitespace.
_rangeExpressionRegex = re.compile(r"\d+(:?\s*-\s*\d+)?(:?\s*,\s*\d+(:?\s*-\s*\d+)?)*")


def parseRangeExpression(s, convertToZeroBased=False):
    """
    Parse a range string expression of the form "1-5,6 & (12 | 100-200)".

    @param s: A C{str} specifiying a range expression, given in the form of
        comma separated numeric ranges or individual indices, interspersed
        with (), |, and &.
    @param convertToZeroBased: If C{True} all indices will have one
        subtracted from them.
    @raise ValueError: If the expression in C{s} cannot be evaluated.
    @return: A C{set} of all C{int}s in the specified expression.
    """
    if not s:
        return set()

    pos = 0
    expr = ""
    for match in _rangeExpressionRegex.finditer(s):
        start, end = match.span()
        if start > pos:
            expr += s[pos:start]
        expr += "{%s}" % ",".join(
            map(str, parseRangeString(match.group(0), convertToZeroBased))
        )
        pos = end

    if pos < len(s):
        expr += s[pos:]

    try:
        return eval(expr)
    except Exception:
        raise ValueError(expr)


def baseCountsToStr(counts):
    """
    Convert base counts to a string.

    @param counts: A C{Counter} instance.
    @return: A C{str} representation of nucleotide counts at an offset.
    """
    return " ".join([("%s:%d" % (base, counts[base])) for base in sorted(counts)])


def nucleotidesToStr(nucleotides, prefix=""):
    """
    Convert offsets and base counts to a string.

    @param nucleotides: A C{defaultdict(Counter)} instance, keyed
        by C{int} offset, with nucleotides keying the Counters.
    @param prefix: A C{str} to put at the start of each line.
    @return: A C{str} representation of the offsets and nucleotide
        counts for each.
    """
    result = []
    for offset in sorted(nucleotides):
        result.append(
            "%s%d: %s" % (prefix, offset, baseCountsToStr(nucleotides[offset]))
        )
    return "\n".join(result)


def countPrint(mesg: str, count: int, len1: int, len2: Optional[int] = None) -> str:
    """
    Format a message followed by an integer count and a percentage (or
    two, if the sequence lengths are unequal).

    @param mesg: a C{str} message.
    @param count: a numeric value.
    @param len1: the C{int} length of sequence 1.
    @param len2: the C{int} length of sequence 2. If C{None}, will
        default to C{len1}.
    @return: A C{str} for printing.
    """

    def percentage(a: int, b: int) -> float:
        """
        What percent of a is b?

        @param a: a numeric value.
        @param b: a numeric value.
        @return: the C{float} percentage.
        """
        return 100.0 * a / b if b else 0.0

    if count == 0:
        return "%s: %d" % (mesg, count)
    else:
        len2 = len2 or len1
        if len1 == len2:
            return "%s: %d/%d (%.2f%%)" % (mesg, count, len1, percentage(count, len1))
        else:
            return (
                "%s: %d/%d (%.2f%%) of sequence 1, "
                "%d/%d (%.2f%%) of sequence 2"
                % (
                    mesg,
                    count,
                    len1,
                    percentage(count, len1),
                    count,
                    len2,
                    percentage(count, len2),
                )
            )


def pct(a: int, b: int) -> str:
    """
    Format a string showing two integers and what percentage the first
    is of the second.

    @param a: An C{int}, the numerator.
    @param b: An C{int}, the denominator.
    """
    assert 0 <= a <= b
    if b:
        return "%d/%d (%.3f%%)" % (a, b, (a / b if b else 0.0) * 100.0)
    else:
        return "0/0 (0.000%)"


@contextmanager
def cd(newdir):
    """
    Trivial context manager for temporarily switching directory.

    @param dir: A C{str} directory to cd to.
    """
    olddir = os.getcwd()
    try:
        os.chdir(newdir)
        yield newdir
    finally:
        os.chdir(olddir)


def take(iterator, n):
    """
    Repeatedly yield lists of C{n} items taken from an iterator.

    @param iterator: Something that can be iterated.
    @param n: An C{int} number of items to yield each time, greater than zero.
    @raise AssertionError: If C{n} < 1.
    @return: A generator that yields lists of length C{n}.
    """
    assert n > 0
    items = []
    iterator = iter(iterator)

    while True:
        try:
            item = next(iterator)
        except StopIteration:
            if items:
                yield items
            break
        else:
            items.append(item)
            if len(items) == n:
                yield items
                items = []


# This function should have a more general name.
def readLabels(fp):
    """
    Read a file of replacement label names.

    @param fp: An open file pointer to read. The file must contain lines with
        a name, a TAB, then a replacement name.
    @return: A C{dict} mapping old names to new names.
    """
    result = {}
    for line in fp:
        oldName, newName = map(str.strip, line.split("\t"))
        result[oldName] = newName
    return result


def matchOffset(leftPaddedReference, leftPaddedQuery):
    """
    At what reference offset does a query begin to match a reference?

    @param leftPaddedReference: A left-padded reference C{str}.
    @param leftPaddedQuery: A left-padded query C{str}.
    @return: An C{int} offset into the (unpadded) reference.
    """
    offset = 0
    for queryChar, referenceChar in zip(leftPaddedQuery, leftPaddedReference):
        if queryChar not in " -":
            break
        offset += referenceChar not in " -"

    return offset


@contextmanager
def openOr(filename, mode="r", defaultFp=None, specialCaseHyphen=True):
    """
    A context manager to either open a file or use a pre-opened default file.

    @param filename: If not C{None}, this is the argument to pass to C{open}
        along with C{mode}. If C{None}, C{fp} is used.
    @param mode: The C{str} file opening mode, used if C{filename} is not
        C{None}.
    @param defaultFp: An open file-like object to yield if C{filename} is
        C{None}.
    @param specialCaseHyphen: If C{True}, treat '-' as a C{None} filename
        and yield the C{defaultFp}.
    """
    if filename is None or filename == "-" and specialCaseHyphen:
        yield defaultFp
    else:
        with open(filename, mode) as fp:
            yield fp

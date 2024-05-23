#!/usr/bin/env python

from math import log10
from typing import Iterator, Optional
from itertools import count

from dark.dimension import dimensionalIterator


def _commonPrefix(a: str, b: str) -> str:
    prefix = []
    for i, j in zip(a, b):
        if i == j:
            prefix.append(i)
        else:
            break

    return "".join(prefix)


def _doInts(
    start: int, end: Optional[int], zeroes: bool, unlimited: bool
) -> Iterator[str]:
    """
    Yield integer ids.
    """
    width = "0" + str(int(log10(end)) + 1) if (zeroes and end is not None) else ""

    if end is None:
        if unlimited:
            # Just produce numbers from start up to (but not including) the
            # next power of ten.
            it = iter(range(start, 10 ** (int(log10(start)) + 1)))
        else:
            # A limit imposed (by the use of maxResults in the call to ids), so it's
            # safe to use itertools.count.
            it = iter(count(start))
    else:
        if start > end:
            start, end = end, start
        it = iter(range(start, end + 1))

    for i in it:
        yield f"{i:{width}d}"


def _doStrings(start: str, end: Optional[str]) -> Iterator[str]:
    if end is not None and len(start) != len(end):
        raise ValueError("Range specifiers must be of equal length.")

    bases = []
    dimensions = []
    odometerStart = []
    parts: list[Optional[str]] = []

    for letter in start:
        if letter.islower():
            first = ord("a")
            bases.append(first)
            dimensions.append(26)
            parts.append(None)
        elif letter.isupper():
            first = ord("A")
            bases.append(first)
            dimensions.append(26)
            parts.append(None)
        elif letter.isdigit():
            first = ord("0")
            bases.append(first)
            dimensions.append(10)
            parts.append(None)
        else:
            # Save this as a constant/unchanging letter seeing as it's not
            # clear how to iterate it.
            parts.append(letter)
            continue

        odometerStart.append(ord(letter) - first)

    for odometer in dimensionalIterator(dimensions, start=odometerStart):
        variable = [f"{chr(bases[i] + c)}" for i, c in enumerate(odometer)]
        resultParts = []
        variableIndex = 0
        for part in parts:
            if part is None:
                resultParts.append(variable[variableIndex])
                variableIndex += 1
            else:
                resultParts.append(part)

        # Make sure we consumed all variable parts from the dimensional iterator.
        assert variableIndex == len(variable)

        result = "".join(resultParts)

        yield result

        if end is not None and result == end:
            break


def ids(
    start: str | int,
    end: Optional[str | int] = None,
    prefix: Optional[str] = None,
    sep: str = " ",
    zeroes: bool = False,
    maxResults: Optional[int] = None,
) -> Iterator[str]:
    """
    Yield ids, starting from C{start}.
    """
    if end is None:
        prefix = prefix or ""
        rest1 = str(start)
        rest2 = None
    else:
        start, end = str(start), str(end)
        if prefix is None:
            prefix = _commonPrefix(start, end)
            rest1 = start[len(prefix) :]
            rest2 = end[len(prefix) :]
        else:
            rest1, rest2 = start, end

    ints = True
    try:
        int1 = int(rest1)
    except ValueError:
        ints = False
    else:
        if rest2 is None:
            int2 = None
        else:
            try:
                int2 = int(rest2)
            except ValueError:
                ints = False

    suffixes = (
        _doInts(int1, int2, zeroes, maxResults is None)
        if ints
        else _doStrings(rest1, rest2)
    )

    resultCount = 0
    for suffix in suffixes:
        if maxResults is not None and resultCount >= maxResults:
            return

        yield f"{prefix}{suffix}"
        resultCount += 1

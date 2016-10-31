def parseBtop(s):
    """
    Parse a btop string.

    @param s: A C{str} btop sequence.
    @raise ValueError: If C{s} is not valid btop.
    @return: A generator that yields integers and 2-tuples of letters, as found
        in C{s}.
    """
    isdigit = str.isdigit
    value = None
    queryLetter = None
    for offset, char in enumerate(s):
        if isdigit(char):
            if queryLetter is not None:
                raise ValueError(
                    'btop string %r has a query letter %r at offset %d with '
                    'no corresponding subject letter' %
                    (s, queryLetter, offset - 1))
            value = int(char) if value is None else value * 10 + int(char)
        else:
            if value is not None:
                yield value
                value = None
                queryLetter = char
            else:
                if queryLetter is None:
                    queryLetter = char
                else:
                    if queryLetter == '-' and char == '-':
                        raise ValueError(
                            'btop string %r has two consecutive gaps at '
                            'offset %d' % (s, offset - 1))
                    elif queryLetter == char:
                        raise ValueError(
                            'btop string %r has two consecutive identical %r '
                            'letters at offset %d' % (s, char, offset - 1))
                    yield (queryLetter, char)
                    queryLetter = None

    if value is not None:
        yield value
    elif queryLetter is not None:
        raise ValueError('btop string %r has a trailing query letter %r with '
                         'no corresponding subject letter' % (s, queryLetter))


def countGaps(s):
    """
    Count the query and subject gaps in a btop string.

    @raises ValueError: If L{parseBtop} finds an error in the btop string C{s}.
    @return: A 2-tuple of C{int}s, with the (query, subject) gaps counts as
        found in C{s}.
    """
    queryGaps = subjectGaps = 0
    for countOrMismatch in parseBtop(s):
        if isinstance(countOrMismatch, tuple):
            queryChar, subjectChar = countOrMismatch
            queryGaps += int(queryChar == '-')
            subjectGaps += int(subjectChar == '-')

    return (queryGaps, subjectGaps)

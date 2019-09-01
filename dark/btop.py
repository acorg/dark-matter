from dark.cigar import CINS, CDEL, CMATCH, CEQUAL, CDIFF


def parseBtop(btopString):
    """
    Parse a BTOP string.

    The format is described at https://www.ncbi.nlm.nih.gov/books/NBK279682/

    @param btopString: A C{str} BTOP sequence.
    @raise ValueError: If C{btopString} is not valid BTOP.
    @return: A generator that yields a series of integers and 2-tuples of
        letters, as found in the BTOP string C{btopString}.
    """
    isdigit = str.isdigit
    value = None
    queryLetter = None
    for offset, char in enumerate(btopString):
        if isdigit(char):
            if queryLetter is not None:
                raise ValueError(
                    'BTOP string %r has a query letter %r at offset %d with '
                    'no corresponding subject letter' %
                    (btopString, queryLetter, offset - 1))
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
                            'BTOP string %r has two consecutive gaps at '
                            'offset %d' % (btopString, offset - 1))
                    elif queryLetter == char:
                        raise ValueError(
                            'BTOP string %r has two consecutive identical %r '
                            'letters at offset %d' %
                            (btopString, char, offset - 1))
                    yield (queryLetter, char)
                    queryLetter = None

    if value is not None:
        yield value
    elif queryLetter is not None:
        raise ValueError(
            'BTOP string %r has a trailing query letter %r with '
            'no corresponding subject letter' % (btopString, queryLetter))


def countGaps(btopString):
    """
    Count the query and subject gaps in a BTOP string.

    @param btopString: A C{str} BTOP sequence.
    @raise ValueError: If L{parseBtop} finds an error in the BTOP string
        C{btopString}.
    @return: A 2-tuple of C{int}s, with the (query, subject) gaps counts as
        found in C{btopString}.
    """
    queryGaps = subjectGaps = 0
    for countOrMismatch in parseBtop(btopString):
        if isinstance(countOrMismatch, tuple):
            queryChar, subjectChar = countOrMismatch
            queryGaps += int(queryChar == '-')
            subjectGaps += int(subjectChar == '-')

    return (queryGaps, subjectGaps)


def btop2cigar(btopString, concise=False, aa=False):
    """
    Convert a BTOP string to a CIGAR string.

    @param btopString: A C{str} BTOP sequence.
    @param concise: If C{True}, use 'M' for matches and mismatches instead
        of the more specific 'X' and '='.
    @param aa: If C{True}, C{btopString} will be interpreted as though it
        refers to amino acids (as in the BTOP string produced by DIAMOND).
        In that case, it is not possible to use the 'precise' CIGAR characters
        because amino acids have multiple codons so we cannot know whether
        an amino acid match is due to an exact nucleotide matches or not.
        Also, the numbers in the BTOP string will be multiplied by 3 since
        they refer to a number of amino acids matching.
    @raise ValueError: If L{parseBtop} finds an error in C{btopString} or
        if C{aa} and C{concise} are both C{True}.
    @return: A C{str} CIGAR string.
    """
    if aa and concise:
        raise ValueError('aa and concise cannot both be True')

    result = []
    thisLength = thisOperation = currentLength = currentOperation = None

    for item in parseBtop(btopString):
        if isinstance(item, int):
            thisLength = item
            thisOperation = CEQUAL if concise else CMATCH
        else:
            thisLength = 1
            query, reference = item
            if query == '-':
                # The query has a gap. That means that in matching the
                # query to the reference a deletion is needed in the
                # reference.
                assert reference != '-'
                thisOperation = CDEL
            elif reference == '-':
                # The reference has a gap. That means that in matching the
                # query to the reference an insertion is needed in the
                # reference.
                thisOperation = CINS
            else:
                # A substitution was needed.
                assert query != reference
                thisOperation = CDIFF if concise else CMATCH

        if thisOperation == currentOperation:
            currentLength += thisLength
        else:
            if currentOperation:
                result.append(
                    '%d%s' %
                    ((3 * currentLength) if aa else currentLength,
                     currentOperation))
            currentLength, currentOperation = thisLength, thisOperation

    # We reached the end of the BTOP string. If there was an operation
    # underway, emit it.  The 'if' here should only be needed to catch the
    # case where btopString was empty.
    assert currentOperation or btopString == ''
    if currentOperation:
        result.append(
            '%d%s' %
            ((3 * currentLength) if aa else currentLength, currentOperation))

    return ''.join(result)

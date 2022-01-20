# From https://samtools.github.io/hts-specs/SAMv1.pdf
CINS, CDEL, CMATCH, CEQUAL, CDIFF = 'IDM=X'


def dna2cigar(s1, s2, concise=False):
    """
    Form a CIGAR string from two equal-length DNA sequences. Case is not
    ignored, so if you want to ignore case you must convert the sequences
    to a consistent case (upper or lower, doesn't matter) before calling.

    @param s1: A C{str} of DNA 'A', 'G', 'C', 'T' characters.
    @param s2: A C{str} of DNA 'A', 'G', 'C', 'T' characters.
    @param concise: If C{True}, use 'M' for matches and mismatches instead
        of the more specific 'X' and '='.
    @raise ValueError: If the two sequences are not of equal length or if the
        sequences are of length zero.
    @return: A C{str} CIGAR string.
    """
    len1 = len(s1)
    if len1 != len(s2):
        raise ValueError('Sequences %r and %r of unequal length (%d != %d).' %
                         (s1, s2, len1, len(s2)))

    if len1 == 0:
        raise ValueError('Two sequences of zero length were passed.')

    if concise:
        return '%d%s' % (len1, CMATCH)

    result = []
    length, operation = 0, None

    for n1, n2 in zip(s1, s2):
        if n1 == n2:
            if operation == CEQUAL:
                length += 1
            else:
                if length:
                    assert operation == CDIFF
                    result.append('%d%s' % (length, CDIFF))
                length, operation = 1, CEQUAL
        else:
            if operation == CDIFF:
                length += 1
            else:
                if length:
                    assert operation == CEQUAL
                    result.append('%d%s' % (length, CEQUAL))
                length, operation = 1, CDIFF

    # Append the final operation.
    result.append('%d%s' % (length, operation))

    return ''.join(result)


def makeCigar(reference, query):
    """
    Make a CIGAR string from an aligned reference and query.

    @param reference: A C{str} reference sequence, possibly padded on
        the left and right with spaces, and possibly with '-' characters
        to indicate that the query has an insertion (relative to the
        reference).
    @param query: A C{str} reference sequence, possibly padded on
        the left with spaces, and possibly with '-' characters to indicate
        that the query has a deletion (relative to the reference).
    @raise ValueError: If the query or reference is empty.
    @return: A C{str} CIGAR string.
    """
    if not reference:
        raise ValueError('Empty reference')
    if not query:
        raise ValueError('Empty query')

    # Pad the reference on the right, if necessary.
    if len(reference) < len(query):
        reference += ' ' * (len(query) - len(reference))

    cigar = []
    softClipLeft = softClipRight = 0
    start = True

    for referenceBase, queryBase in zip(reference, query):
        if referenceBase == ' ':
            if queryBase == ' ':
                continue
            else:
                if start:
                    softClipLeft += 1
                else:
                    softClipRight += 1
        else:
            start = False
            if queryBase == ' ':
                continue
            elif referenceBase == '-':
                # Insertion to the reference.
                assert queryBase != '-'
                cigar.append(CINS)
            elif queryBase == '-':
                # Deletion from the reference.
                assert referenceBase != '-'
                cigar.append(CDEL)
            else:
                cigar.append(CMATCH)

    # Replace strings of identical operations with a count and the operation.
    lastOp = None
    count = 0
    middle = []
    for op in cigar:
        if op != lastOp:
            if count:
                middle.append(f'{count}{lastOp}')
            count = 1
            lastOp = op
        else:
            count += 1

    if count:
        middle.append(f'{count}{lastOp}')

    return ((f'{softClipLeft}S' if softClipLeft else '') +
            ''.join(middle) +
            (f'{softClipRight}S' if softClipRight else ''))

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

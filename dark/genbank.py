def splitRange(rangeStr):
    """
    Split a GenBank range string.

    @param rangeStr: A C{str} indicating the (1-based) genome nucleotide
        offsets covered by a protein. This can have the following forms:

            'A..B' where A and B are C{int} offsets, with A <= B.
            'MOD(A..B,C..D,...)' where A, B, etc. are pairs of non-decreasing
                integers, as above and MOD is either 'join' or 'complement'.
            'complement(join(A..B,C..D,...))' with A, B, etc. as above. When
                a range is both a join and a complement, the two modifiers
                only occur in this order.

        Note that this is the string as found in a GenBank flat file, not the
        string that is returned from BioPython SeqIO.parse when it parses such
        a file. For that use C{splitBioPythonRange}, below.

    @raise ValueError: If C{rangeStr} cannot be correctly parsed.
    @return: A 2-tuple consisting of a C{tuple} and a C{bool}. The C{tuple}
            consists of 2-C{tuples} of C{int}s (following the Python string
            offset convention) corresponding to the A, B, etc., above. The
            C{bool} indicates whether 'complement' was found in the range
            string. Example arguments and results are:

                '3..5'                          -> (((3, 6)), False)
                'join(3..5,70..77)'             -> (((3, 6), (70, 78)), False)
                'complement(join(3..5,70..77))' -> (((3, 6), (70, 78)), True)
    """

    if rangeStr.startswith('complement(') and rangeStr[-1] == ')':
        complement = True
        inner = rangeStr[11:-1]
    else:
        complement = False
        inner = rangeStr

    if inner.startswith('join(') and inner[-1] == ')':
        join = True
        inner = inner[5:-1]
    else:
        join = False

    ranges = []
    subRanges = inner.split(',')

    if len(subRanges) == 1 and join:
        raise ValueError('Could not parse GenBank range string "%s". '
                         'join() can only be used with multiple ranges.' %
                         rangeStr)
    elif len(subRanges) > 1 and not join:
        raise ValueError('Could not parse GenBank range string "%s". '
                         'Multiple ranges must be wrapped in join().' %
                         rangeStr)

    for subRange in subRanges:
        try:
            start, stop = map(int, subRange.split('..'))
        except ValueError as e:
            raise ValueError('Could not parse GenBank range string "%s". '
                             'Original parsing ValueError was "%s".' %
                             (rangeStr, e))
        else:
            if start > stop:
                raise ValueError('Could not parse GenBank range string "%s". '
                                 'Offset values (%d, %d) cannot decrease.' %
                                 (rangeStr, start, stop))
            # Convert to Python offsets.
            ranges.append((start - 1, stop))

    return (tuple(ranges), complement)


def splitBioPythonRange(rangeStr):
    """
    Split a GenBank range string as converted by BioPython.

    @param rangeStr: A C{str} indicating the (0-based) genome nucleotide
        offsets covered by a protein. This can have the following example
        forms:

        [9462:10137](+)
        [11969:12575](-)
        join{[126386:126881](-), [125941:126232](+)}
        join{[153005:153821](+), [82030:82618](-), [80480:81305](-)}

        Note that this is the string as returned from BioPython SeqIO.parse
        when it parses a GenBank flat file, not the string that is in such a
        file. For that use C{splitRange}, above.

    @raise ValueError: If C{rangeStr} cannot be correctly parsed.
    @return: A C{list} whose elements are 3-C{tuple}s holding the C{int} start
        and stop offsets and a C{bool} that indicates whether the offset
        corresponds to the complement strand (i.e., this will be C{False} when
        there is a '(+)' in a range and C{True} when there is a '(-)'. Example
        arguments and their return values:

            '[9462:10137](+)' => ((9462, 10137, False))

            '[11969:12575](-)' => ((11969, 12575, True))

            'join{[126386:126881](-), [125941:126232](+)}'  =>
                ((126386, 126881, True), (125941, 126232, False))
    """

    if rangeStr.startswith('join{') and rangeStr[-1] == '}':
        join = True
        inner = rangeStr[5:-1]
    else:
        join = False
        inner = rangeStr

    ranges = []
    subRanges = inner.split(', ')

    if len(subRanges) == 1 and join:
        raise ValueError('Could not parse GenBank range string "%s". '
                         'join{} can only be used with multiple ranges.' %
                         rangeStr)
    elif len(subRanges) > 1 and not join:
        raise ValueError('Could not parse GenBank range string "%s". '
                         'Multiple ranges must be wrapped in join{}.' %
                         rangeStr)

    for subRange in subRanges:
        if subRange.endswith('](+)'):
            complement = False
        elif subRange.endswith('](-)'):
            complement = True
        else:
            raise ValueError('Could not parse GenBank range string "%s". '
                             'Range "%s" does not end with ](+) or ](-).' %
                             (rangeStr, subRange))

        if not subRange.startswith('['):
            raise ValueError('Could not parse GenBank range string "%s". '
                             'Range "%s" does not start with "[".' %
                             (rangeStr, subRange))
        try:
            start, stop = subRange[1:-4].split(':')
        except ValueError as e:
            raise ValueError('Could not parse GenBank range string "%s". '
                             'Original parsing ValueError was "%s".' %
                             (rangeStr, e))
        else:
            if start.startswith('<'):
                start = start[1:]
            if stop.startswith('>'):
                stop = stop[1:]
            start, stop = map(int, (start, stop))
            if start > stop:
                raise ValueError('Could not parse GenBank range string "%s". '
                                 'Offset values (%d, %d) cannot decrease.' %
                                 (rangeStr, start, stop))
            ranges.append((start, stop, complement))

    return ranges

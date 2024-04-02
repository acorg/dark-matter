from pysam import CINS, CSOFT_CLIP, CHARD_CLIP  # type: ignore

from dark.sam import CONSUMES_REFERENCE


# From https://samtools.github.io/hts-specs/SAMv1.pdf
(CINS_STR, CDEL_STR, CMATCH_STR, CEQUAL_STR, CDIFF_STR, CHARD_CLIP_STR) = tuple(
    "IDM=XH"
)


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
        raise ValueError(
            "Sequences %r and %r of unequal length (%d != %d)."
            % (s1, s2, len1, len(s2))
        )

    if len1 == 0:
        raise ValueError("Two sequences of zero length were passed.")

    if concise:
        return "%d%s" % (len1, CMATCH_STR)

    result = []
    length, operation = 0, None

    for n1, n2 in zip(s1, s2):
        if n1 == n2:
            if operation == CEQUAL_STR:
                length += 1
            else:
                if length:
                    assert operation == CDIFF_STR
                    result.append("%d%s" % (length, CDIFF_STR))
                length, operation = 1, CEQUAL_STR
        else:
            if operation == CDIFF_STR:
                length += 1
            else:
                if length:
                    assert operation == CEQUAL_STR
                    result.append("%d%s" % (length, CEQUAL_STR))
                length, operation = 1, CDIFF_STR

    # Append the final operation.
    result.append("%d%s" % (length, operation))

    return "".join(result)


def makeCigar(reference, query, noEdgeInsertions=True):
    """
    Make a CIGAR string from an aligned reference and query.

    @param reference: A C{str} reference sequence, possibly padded on
        the left and right with spaces, and possibly with '-' characters
        to indicate that the query has an insertion (relative to the
        reference).
    @param query: A C{str} reference sequence, possibly padded on
        the left with spaces, and possibly with '-' characters to indicate
        that the query has a deletion (relative to the reference).
    @param noEdgeInsertions: If {True}, convert a leading and trailing
        insertions in the result into an equivalent number of soft clips.
        This is useful when making a CIGAR string that would otherwise have
        insertions due to the reference containing '-' characters.
    @raise ValueError: If the query or reference is empty.
    @return: A C{str} CIGAR string.
    """
    if not reference:
        raise ValueError("Empty reference")
    if not query:
        raise ValueError("Empty query")

    # Pad the reference on the right, if necessary.
    if len(reference) < len(query):
        reference += " " * (len(query) - len(reference))

    cigar = []
    softClipLeft = softClipRight = 0
    start = True

    for referenceBase, queryBase in zip(reference, query):
        if referenceBase == " ":
            if queryBase == " ":
                continue
            else:
                if start:
                    softClipLeft += 1
                else:
                    softClipRight += 1
        else:
            start = False
            if queryBase == " ":
                continue
            elif referenceBase == "-":
                # Insertion to the reference.
                assert queryBase != "-"
                cigar.append(CINS_STR)
            elif queryBase == "-":
                # Deletion from the reference.
                assert referenceBase != "-"
                cigar.append(CDEL_STR)
            else:
                cigar.append(CMATCH_STR)

    # Replace strings of identical operations with a count and the operation.
    lastOp = None
    count = 0
    middle = []
    for op in cigar:
        if op != lastOp:
            if count:
                middle.append(f"{count}{lastOp}")
            count = 1
            lastOp = op
        else:
            count += 1

    if count:
        middle.append(f"{count}{lastOp}")

    if noEdgeInsertions:
        if middle:
            if middle[0].endswith(CINS_STR):
                softClipLeft += int(middle[0][:-1])
                middle.pop(0)

        if middle:
            if middle[-1].endswith(CINS_STR):
                softClipRight += int(middle[-1][:-1])
                middle.pop()

    return (
        (f"{softClipLeft}S" if softClipLeft else "")
        + "".join(middle)
        + (f"{softClipRight}S" if softClipRight else "")
    )


def cigarTuplesToOperations(tuples, includeHardClip=True):
    """
    Produce a sequence of CIGAR operations given a cigar tuples list
    from a pysam AlignedRead instance.

    See https://pysam.readthedocs.io/en/latest/
        api.html#pysam.AlignedSegment.cigartuples

    @param tuples: A C{list} of 2-tuples, each containing an operation
        and an C{int} count.
    @param includeHardClip: If C{True}, yield hard clipping operations,
        else elide them.
    """
    for operation, count in tuples:
        if operation != CHARD_CLIP or includeHardClip:
            for _ in range(count):
                yield operation


def softClippedOffset(offset, pairs, cigarOperations):
    """
    Get the offset in the reference of a base that has been soft-clipped
    in a query.

    @param offset: The offset (in C{pairs}) of the soft-clipped base (which
       has a to-be-determined offset in the reference).
    @param pairs: A C{list} of C{int} (queryOffset, referenceOffset) pairs
       for a pysam C{pysam.AlignedSegment}, as returned by
       C{read.get_aligned_pairs} where the read is returned by the pysam
       C{fetch} method.
    @param cigarOperations: A C{list} of pysam CIGAR operations.
    @return: The C{int} reference offset of the soft-clipped base.
    """
    assert cigarOperations[offset] == CSOFT_CLIP
    assert pairs[offset][1] is None

    # Look back.
    count = 0
    for pair, cigarOperation in zip(pairs[offset::-1], cigarOperations[offset::-1]):
        if cigarOperation == CSOFT_CLIP:
            count += 1
        else:
            _, referenceOffset = pair
            # If we have a reference offset, then we've found our
            # answer. Otherwise, this has to be a hard-clipped base.
            if referenceOffset is None:
                assert cigarOperation == CHARD_CLIP
                # We could break but won't, as explained above.
            else:
                return referenceOffset + count

    # Look ahead.
    count = 0
    for pair, cigarOperation in zip(pairs[offset:], cigarOperations[offset:]):
        if cigarOperation == CSOFT_CLIP:
            count += 1
        else:
            _, referenceOffset = pair
            # If we have a reference offset, then we've found our
            # answer. Otherwise, this has to be a hard-clipped base.
            if referenceOffset is None:
                assert cigarOperation == CHARD_CLIP
                # We could break here, but if we keep going we can test
                # that the whole rest of the CIGAR is hard clips and will
                # fail if not. That's not needed in a world of perfect code.
            else:
                return referenceOffset - count

    # This should be impossible.
    raise ValueError(
        "Soft-clipped base with no following or preceding non-hard-clipped bases."
    )


def insertionOffset(offset, pairs, cigarOperations):
    """
    Get the insertion offset in the reference for a base that is part of an
    insertion.

    @param offset: The offset (in C{pairs}) of the soft-clipped base (which
       has a to-be-determined offset in the reference).
    @param pairs: A C{list} of C{int} (queryOffset, referenceOffset) pairs
       for a pysam C{pysam.AlignedSegment}, as returned by
       C{read.get_aligned_pairs} where the read is returned by the pysam
       C{fetch} method.
    @param cigarOperations: A C{list} of pysam CIGAR operations.
    @return: A 2-tuple containing a C{bool} and the C{int} reference offset of
        the base after the insertion. The C{bool} indicates whether we needed
        to look back in the sequence to find a reference base.
    """
    assert cigarOperations[offset] == CINS
    assert pairs[offset][1] is None

    # Look back.
    for pair, cigarOperation in zip(pairs[offset::-1], cigarOperations[offset::-1]):
        if cigarOperation != CINS:
            _, referenceOffset = pair
            # If we have a reference offset, then we've found our
            # answer. Otherwise, check as above.
            if referenceOffset is None:
                assert cigarOperation not in CONSUMES_REFERENCE
            else:
                return True, referenceOffset + 1

    # Look ahead.
    for pair, cigarOperation in zip(pairs[offset:], cigarOperations[offset:]):
        if cigarOperation != CINS:
            _, referenceOffset = pair
            # If we have a reference offset, then we've found our
            # answer. Otherwise, check that this is something that does not
            # consume reference bases.
            if referenceOffset is None:
                assert cigarOperation not in CONSUMES_REFERENCE
            else:
                return False, referenceOffset

    # This should be impossible.
    raise ValueError("Inserted base with no following or preceding reference bases.")

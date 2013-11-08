def normalizeHSP(hsp, queryLen):
    """
    Examime the sense of an HSP and return information about where the
    query and the alignment (match) begin and end.  BLAST always returns
    query start and stop values that are increasing, but the reported
    subject match may be reversed (start > stop).  Return a dict with keys
    that allow the query and the alignment to be displayed relative to the
    subject orientation (i.e., with start < stop for both the query and the
    match).

    In the returned object, all indices are suitable for Python string
    slicing etc.  We must be careful to convert from the 1-based offsets
    found in BLAST output properly.  The values in the returned dictionary
    are ALL indices into the subject string.

    hsp.frame is a (query, subject) 2-tuple, with both values coming from
    the set {-3, -2, -1, 1, 2, 3}. The query value indicates negative or
    positive sense (negative or positive sign). The subject value is the
    starting nucleotide match offset modulo 3, plus one (i.e., it tells us
    which of the 3 possible reading frames is used in the match. It is
    redundant because that information can also be obtained from the mod 3
    value of the match offset). In pure nucleotide matching, the

    NOTE: the returned queryStart value may be negative.  The subject
    sequence is considered to start at offset 0.  So if the query string
    has sufficient additional nucleotides before the start of the alignment
    match, it may protrude to the left of the subject (i.e., have x < 0).

    hsp: a HSP from a BLAST record.  All passed hsp offsets are 1-based.
    queryLen: the length of the query sequence.

    """

    queryPositive = hsp.frame[0] > 0
    subjectPositive = hsp.frame[1] > 0

    query_start = hsp.query_start
    query_end = hsp.query_end
    sbjct_start = hsp.sbjct_start
    sbjct_end = hsp.sbjct_end

    # Query indices are always ascending.
    assert query_start <= query_end

    # Subject indices are only descending when the subject sense is negative.
    if sbjct_start > sbjct_end:
        assert not subjectPositive
        # Make subject indices ascending.
        sbjct_start, sbjct_end = sbjct_end, sbjct_start

    # Sanity check that the length of the matches in the subject and query
    # are identical, taking into account gaps in either (indicated by '-'
    # characters in the match sequences, as returned by BLAST).
    subjectLength = sbjct_end - sbjct_start + 1 + hsp.sbjct.count('-')
    queryLength = query_end - query_start + 1 + hsp.query.count('-')
    assert subjectLength == queryLength, (
        'Subject match length (%d) != Query match length (%d)' %
        (subjectLength, queryLength))

    # TODO: check the mod 3 value of the start offsets.

    # Set subject indices to be zero-based.
    subjectStart = sbjct_start - 1
    subjectEnd = sbjct_end

    # Calculate query indices. These are indices relative to the subject!

    # unmatchedQueryLeft is the number of query bases that will be sticking
    # out to the left of the start of the subject in our plots.
    if queryPositive:
        unmatchedQueryLeft = query_start - 1
    else:
        unmatchedQueryLeft = queryLen - query_end

    # Set the query offsets based on the direction the match with the
    # subject takes.
    if subjectPositive:
        queryStart = subjectStart - unmatchedQueryLeft
        queryEnd = queryStart + queryLen
    else:
        queryEnd = subjectEnd + unmatchedQueryLeft
        queryStart = queryEnd - queryLen

    # Final sanity check.
    assert queryStart <= subjectStart
    assert queryEnd >= subjectEnd

    return {
        'subjectEnd': subjectEnd,
        'subjectStart': subjectStart,
        'queryEnd': queryEnd,
        'queryStart': queryStart,
    }

from math import exp


def printHSP(hsp, indent=''):
    for attr in ['align_length', 'bits', 'expect', 'frame', 'gaps',
                 'identities', 'num_alignments', 'positives', 'query_end',
                 'query_start', 'sbjct', 'match', 'query', 'sbjct_end',
                 'sbjct_start', 'score', 'strand']:
        print '%s%s: %s' % (indent, attr, getattr(hsp, attr))
    print '%sp: %.10f' % (indent, 1.0 - exp(-1.0 * hsp.expect))


def normalizeHSP(hsp, queryLen):
    """
    Examime the sense of an HSP and return information about where the
    query and the alignment (match) begin and end.  Return a dict with keys
    that allow the query and the alignment to be displayed relative to the
    subject orientation (i.e., with start < stop for both the query and the
    match). The returned query indices are offsets into the subject. I.e.,
    they indicate where on the subject the query lies.

    In the returned object, all indices are suitable for Python string
    slicing etc.  We must be careful to convert from the 1-based offsets
    found in BLAST output properly.

    hsp.frame is a (query, subject) 2-tuple, with both values coming from
    {-3, -2, -1, 1, 2, 3}. The sign indicates negative or positive sense
    (i.e., the direction of reading through the query or subject to get the
    alignment). The value is the nucleotide match offset modulo 3, plus one
    (i.e., it tells us which of the 3 possible reading frames is used in
    the match). The value is redundant because that information could also
    be obtained from the mod 3 value of the match offset.

    NOTE: the returned queryStart value may be negative.  We consider the
    subject sequence to start at offset 0.  So if the query string has
    sufficient additional nucleotides before the start of the alignment
    match, it may protrude to the left of the subject. Similarly, the
    returned queryEnd can be greater than the subjectEnd.

    hsp: a HSP from a BLAST record.  All passed hsp offsets are 1-based.
    queryLen: the length of the query sequence.
    """

    queryPositive = hsp.frame[0] > 0
    subjectPositive = hsp.frame[1] > 0

    # The following variable names, with underscores, match the names of
    # attributes BioPython uses and the values (1-based) match those
    # reported by BLAST.
    query_start = hsp.query_start
    query_end = hsp.query_end
    sbjct_start = hsp.sbjct_start
    sbjct_end = hsp.sbjct_end

    # When the query is positive, BLASTN and TBLASTX give query offsets
    # ascending.
    #
    # TBLASTX reports negative query sense with indices ascending.
    # BLASTN does not report negative query sense.
    #
    # In all cases the query offsets should be ascending.
    assert query_start <= query_end, (
        'Assertion "query_start <= query_end" failed. Query positive is %s. '
        'query_start = %d, query_end = %d' %
        (queryPositive, query_start, query_end))

    if subjectPositive:
        # Make sure indices are ascending.
        assert sbjct_start <= sbjct_end
    else:
        # Subject is negative. Its indices will be ascending for TBLASTX
        # output but descending for BLASTN :-( Make sure we have them
        # ascending.
        if sbjct_start > sbjct_end:
            sbjct_start, sbjct_end = sbjct_end, sbjct_start

    # Sanity check that the length of the matches in the subject and query
    # are identical, taking into account gaps in either (indicated by '-'
    # characters in the match sequences, as returned by BLAST).
    subjectGaps = hsp.sbjct.count('-')
    queryGaps = hsp.query.count('-')
    subjectLengthWithGaps = sbjct_end - sbjct_start + 1 + subjectGaps
    queryLengthWithGaps = query_end - query_start + 1 + queryGaps
    assert subjectLengthWithGaps == queryLengthWithGaps, (
        'Including gaps, subject match length (%d) != Query match length (%d)'
        % (subjectLengthWithGaps, queryLengthWithGaps))

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
        queryEnd = queryStart + queryLen + queryGaps
    else:
        queryEnd = subjectEnd + unmatchedQueryLeft

        queryStart = queryEnd - queryLen - queryGaps

    # TODO: remove all this once we're sure we have things right.
    #
    # def debugPrint():
    #     print 'frame', hsp.frame
    #     print 'queryLen is', queryLen
    #     print 'queryPositive: %s' % queryPositive
    #     print 'subjectPositive: %s' % subjectPositive
    #     print 'unmatchedQueryLeft: %s' % unmatchedQueryLeft
    #     print 'subjectGaps: %s' % subjectGaps
    #     print 'queryGaps: %s' % queryGaps
    #     print 'query_start: %s' % query_start
    #     print 'query_end: %s' % query_end
    #     print 'sbjct_start: %s' % sbjct_start
    #     print 'sbjct_end: %s' % sbjct_end
    #     print 'queryStart: %s' % queryStart
    #     print 'queryEnd: %s' % queryEnd
    #     print 'subjectStart: %s' % subjectStart
    #     print 'subjectEnd: %s' % subjectEnd
    #     print 'query  : %s' % hsp.query
    #     print 'subject: %s' % hsp.sbjct
    #
    # if queryStart > subjectStart:
    #     print 'Oops: queryStart > subjectStart'
    #     debugPrint()
    # if queryEnd < subjectEnd:
    #     print 'Oops: queryEnd < subjectEnd'
    #     debugPrint()

    # Final sanity checks.
    #assert queryStart <= subjectStart
    #assert queryEnd >= subjectEnd

    return {
        'queryEnd': queryEnd,
        'queryStart': queryStart,
        'subjectEnd': subjectEnd,
        'subjectStart': subjectStart,
    }

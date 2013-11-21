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

    hsp.frame is a (query, subject) 2-tuple, with the query value coming
    from {1, 2, 3} and the subject from {-3, -2, -1, 1, 2, 3}. The sign
    indicates negative or positive sense (negative or positive sign). The
    subject value is the starting nucleotide match offset modulo 3, plus
    one (i.e., it tells us which of the 3 possible reading frames is used
    in the match. It is redundant because that information can also be
    obtained from the mod 3 value of the match offset.

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

    # As far as we know (blastn and tblastx), the query frame is always
    # positive and the query indices are always ascending.
    assert queryPositive and query_start <= query_end, (
        'Assertion "queryPositive and query_start <= query_end" failed. '
        'queryPositive = %s, query_start = %d, query_end = %d' % (
            queryPositive, query_start, query_end))

    if subjectPositive:
        # Make sure indices are ascending.
        assert sbjct_start <= sbjct_end
    else:
        # Subject is negative. Its indices will be ascending for tblastx output
        # but descending for blastn :-(  Make sure we have them ascending.
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
    assert queryStart <= subjectStart
    assert queryEnd >= subjectEnd

    return {
        'queryEnd': queryEnd,
        'queryStart': queryStart,
        'subjectEnd': subjectEnd,
        'subjectStart': subjectStart,
    }

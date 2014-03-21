from math import exp


def printHSP(hsp, indent=''):
    for attr in ['align_length', 'bits', 'expect', 'frame', 'gaps',
                 'identities', 'num_alignments', 'positives', 'query_end',
                 'query_start', 'sbjct', 'match', 'query', 'sbjct_end',
                 'sbjct_start', 'score', 'strand']:
        print '%s%s: %s' % (indent, attr, getattr(hsp, attr))
    print '%sp: %.10f' % (indent, 1.0 - exp(-1.0 * hsp.expect))


def normalizeHSP(hsp, queryLen, blastApplication):
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

    @param hsp: a HSP from a BLAST record.  All passed hsp offsets are 1-based.
    @param queryLen: the length of the query sequence.
    @param blastApplication: The C{str} command line program that was
        run (e.g., 'blastn', 'blastx').
    """

    def debugPrint(locals, msg=None):
        """
        Print debugging information showing the local variables from
        a call to normalizeHSP and then raise an C{AssertionError}.

        @param locals: A C{dict} of local variables.
        @param msg: A C{str} message to raise C{AssertionError} with.
        """
        print 'normalizeHSP error:'
        print '  queryLen: %d' % queryLen
        for var in sorted(locals.keys()):
            if var in ('debugPrint', 'hsp'):
                continue
            print '  %s: %s' % (var, locals[var])
        print '  Original HSP:'
        printHSP(hsp, '    ')
        if msg:
            raise AssertionError(msg)
        else:
            raise AssertionError()

    queryPositive = hsp.frame[0] > 0
    subjectPositive = hsp.frame[1] > 0

    # The following variable names with underscores match the names of
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
    if query_start > query_end:
        debugPrint(locals(),
                   'Assertion "query_start <= query_end" failed. Query '
                   'positive is %s. query_start = %d, query_end = %d' %
                   (queryPositive, query_start, query_end))

    if subjectPositive:
        # Make sure indices are ascending.
        if sbjct_start > sbjct_end:
            debugPrint(locals())
    else:
        # Subject is negative. Its indices will be ascending for TBLASTX
        # output but descending for BLASTN :-( Make sure we have them
        # ascending.
        if sbjct_start > sbjct_end:
            sbjct_start, sbjct_end = sbjct_end, sbjct_start

    # Now that we have asserted what we can about the original HSP values
    # and gotten them into ascending order, make some sane 0-based offsets.
    queryStart = query_start - 1
    queryEnd = query_end
    subjectStart = sbjct_start - 1
    subjectEnd = sbjct_end

    if blastApplication == 'blastx':
        # In Blastx output, subject offsets are based on protein sequence
        # length but queries (and the reported offsets) are nucleotide.
        # Convert the query offsets to protein because we will plot against
        # the subject (protein).
        #
        # Note that queryStart and queryEnd may not be 0 mod 3. They are
        # offsets into the query string giving the position of the AA,
        # which depends on the translation frame.
        queryStart = int(queryStart / 3)
        queryEnd = int(queryEnd / 3)

    # No operations on original 1-based HSP variables (with underscores)
    # should appear beyond this point.

    subjectLength = subjectEnd - subjectStart
    queryLength = queryEnd - queryStart

    subjectGaps = hsp.sbjct.count('-')
    queryGaps = hsp.query.count('-')

    # Sanity check that the length of the matches in the subject and query
    # are identical, taking into account gaps in either (indicated by '-'
    # characters in the match sequences, as returned by BLAST).
    subjectLengthWithGaps = subjectLength + subjectGaps
    queryLengthWithGaps = queryLength + queryGaps
    if subjectLengthWithGaps != queryLengthWithGaps:
        debugPrint(locals(),
                   'Including gaps, subject match length (%d) != Query match '
                   'length (%d)' % (subjectLengthWithGaps,
                                    queryLengthWithGaps))

    # TODO: check the mod 3 value of the start offsets.

    # Calculate query indices. These are indices relative to the subject!

    # unmatchedQueryLeft is the number of query bases that will be sticking
    # out to the left of the start of the subject in our plots.
    if queryPositive:
        unmatchedQueryLeft = queryStart
    else:
        unmatchedQueryLeft = queryLen - queryEnd

    # Set the query offsets based on the direction the match with the
    # subject takes.
    if subjectPositive:
        queryStart = subjectStart - unmatchedQueryLeft
        queryEnd = queryStart + queryLen + queryGaps
    else:
        queryEnd = subjectEnd + unmatchedQueryLeft

        queryStart = queryEnd - queryLen - queryGaps

    # Final sanity checks.
    if queryStart > subjectStart:
        debugPrint(locals(), 'queryStart > subjectStart')
    if queryEnd < subjectEnd:
        debugPrint(locals(), 'queryEnd < subjectEnd')

    return {
        'queryEnd': queryEnd,
        'queryStart': queryStart,
        'subjectEnd': subjectEnd,
        'subjectStart': subjectStart,
    }

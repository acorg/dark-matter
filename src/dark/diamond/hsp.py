import sys

from dark.btop import countGaps


def _debugPrint(hsp, queryLen, localDict, msg=""):
    """
    Print debugging information showing the local variables used during
    a call to normalizeHSP and the hsp and then raise an C{AssertionError}.

    @param hsp: The HSP C{dict} passed to normalizeHSP.
    @param queryLen: the length of the query sequence.
    @param localDict: A C{dict} of local variables (as produced by locals()).
    @param msg: A C{str} message to raise C{AssertionError} with.
    @raise AssertionError: unconditionally.
    """
    print("normalizeHSP error:", file=sys.stderr)
    print("  queryLen: %d" % queryLen, file=sys.stderr)

    print("  Original HSP:", file=sys.stderr)
    for attr in [
        "bits",
        "btop",
        "expect",
        "frame",
        "query_end",
        "query_start",
        "sbjct",
        "query",
        "sbjct_end",
        "sbjct_start",
    ]:
        print("    %s: %r" % (attr, hsp[attr]), file=sys.stderr)

    print("  Local variables:", file=sys.stderr)
    for var in sorted(localDict):
        if var != "hsp":
            print("    %s: %s" % (var, localDict[var]), file=sys.stderr)

    raise AssertionError(msg)


def _sanityCheck(
    subjectStart,
    subjectEnd,
    queryStart,
    queryEnd,
    queryStartInSubject,
    queryEndInSubject,
    hsp,
    queryLen,
    subjectGaps,
    queryGaps,
    localDict,
):
    """
    Perform some sanity checks on an HSP. Call _debugPrint on any error.

    @param subjectStart: The 0-based C{int} start offset of the match in the
        subject.
    @param subjectEnd: The 0-based C{int} end offset of the match in the
        subject.
    @param queryStart: The 0-based C{int} start offset of the match in the
        query.
    @param queryEnd: The 0-based C{int} end offset of the match in the query.
    @param queryStartInSubject: The 0-based C{int} offset of where the query
        starts in the subject.
    @param queryEndInSubject: The 0-based C{int} offset of where the query
        ends in the subject.
    @param hsp: The HSP C{dict} passed to normalizeHSP.
    @param queryLen: the C{int} length of the query sequence.
    @param subjectGaps: the C{int} number of gaps in the subject.
    @param queryGaps: the C{int} number of gaps in the query.
    @param localDict: A C{dict} of local variables from our caller (as
        produced by locals()).
    """
    # Subject indices must always be ascending.
    if subjectStart >= subjectEnd:
        _debugPrint(hsp, queryLen, localDict, "subjectStart >= subjectEnd")

    subjectMatchLength = subjectEnd - subjectStart
    queryMatchLength = queryEnd - queryStart

    # Sanity check that the length of the matches in the subject and query
    # are identical, taking into account gaps in both.
    subjectMatchLengthWithGaps = subjectMatchLength + subjectGaps
    queryMatchLengthWithGaps = queryMatchLength + queryGaps
    if subjectMatchLengthWithGaps != queryMatchLengthWithGaps:
        _debugPrint(
            hsp,
            queryLen,
            localDict,
            "Including gaps, subject match length (%d) != Query match "
            "length (%d)" % (subjectMatchLengthWithGaps, queryMatchLengthWithGaps),
        )

    if queryStartInSubject > subjectStart:
        _debugPrint(
            hsp,
            queryLen,
            localDict,
            "queryStartInSubject (%d) > subjectStart (%d)"
            % (queryStartInSubject, subjectStart),
        )
    if queryEndInSubject < subjectEnd:
        _debugPrint(
            hsp,
            queryLen,
            localDict,
            "queryEndInSubject (%d) < subjectEnd (%d)"
            % (queryEndInSubject, subjectEnd),
        )


def normalizeHSP(hsp, queryLen, diamondTask):
    """
    Examine an HSP and return information about where the query and subject
    match begins and ends.  Return a dict with keys that allow the query to
    be displayed against the subject. The returned readStartInSubject and
    readEndInSubject indices are offsets into the subject. I.e., they
    indicate where in the subject the query falls.

    In the returned object, all indices are suitable for Python string
    slicing etc.  We must be careful to convert from the 1-based offsets
    found in DIAMOND output properly.

    hsp['frame'] is a value from {-3, -2, -1, 1, 2, 3}. The sign indicates
    negative or positive sense (i.e., the direction of reading through the
    query to get the alignment). The frame value is the nucleotide match offset
    modulo 3, plus one (i.e., it tells us which of the 3 possible query reading
    frames was used in the match).

    NOTE: the returned readStartInSubject value may be negative.  We consider
    the subject sequence to start at offset 0.  So if the query string has
    sufficient additional nucleotides before the start of the alignment
    match, it may protrude to the left of the subject. Similarly, the returned
    readEndInSubject can be greater than the subjectEnd.

    @param hsp: an HSP in the form of a C{dict}, built from a DIAMOND record.
        All passed offsets are 1-based.
    @param queryLen: the length of the query sequence.
    @param diamondTask: The C{str} command-line matching algorithm that was
        run (either 'blastx' or 'blastp').
    @return: A C{dict} with C{str} keys and C{int} offset values. Keys are
            readStart
            readEnd
            readStartInSubject
            readEndInSubject
            subjectStart
            subjectEnd
        The returned offset values are all zero-based.
    """

    queryGaps, subjectGaps = countGaps(hsp["btop"])

    # Make some variables using Python's standard string indexing (start
    # offset included, end offset not). No calculations in this function
    # are done with the original 1-based HSP variables.
    queryStart = hsp["query_start"] - 1
    queryEnd = hsp["query_end"]
    subjectStart = hsp["sbjct_start"] - 1
    subjectEnd = hsp["sbjct_end"]

    queryReversed = hsp["frame"] < 0

    # Query offsets must be ascending, unless we're looking at blastx output
    # and the query was reversed for the match.
    if queryStart >= queryEnd:
        if diamondTask == "blastx" and queryReversed:
            # Compute new query start and end indices, based on their
            # distance from the end of the string.
            #
            # Above we took one off the start index, so we need to undo
            # that (because the start is actually the end). We didn't take
            # one off the end index, and need to do that now (because the
            # end is actually the start).
            queryStart = queryLen - (queryStart + 1)
            queryEnd = queryLen - (queryEnd - 1)
        else:
            _debugPrint(hsp, queryLen, locals(), "queryStart >= queryEnd")

    if diamondTask == "blastx":
        # In DIAMOND blastx output, subject offsets are based on protein
        # sequence length but queries (and the reported offsets) are
        # nucleotide.  Convert the query offsets to protein because we will
        # plot against the subject (protein).
        #
        # Convert queryLen and the query nucleotide start and end offsets
        # to be valid for the query after translation to AAs. When
        # translating, DIAMOND may ignore some nucleotides at the start
        # and/or the end of the original DNA query. At the start this is
        # due to the frame in use, and at the end it is due to always using
        # three nucleotides at a time to form codons.
        #
        # So, for example, a query of 6 nucleotides that is translated in
        # frame 2 (i.e., the translation starts from the second nucleotide)
        # will have length 1 as an AA sequence. The first nucleotide is
        # ignored due to the frame and the last two due to there not being
        # enough final nucleotides to make another codon.
        #
        # In the following, the subtraction accounts for the first form of
        # loss and the integer division for the second.
        initiallyIgnored = abs(hsp["frame"]) - 1
        queryLen = (queryLen - initiallyIgnored) // 3
        queryStart = (queryStart - initiallyIgnored) // 3
        queryEnd = (queryEnd - initiallyIgnored) // 3

    # unmatchedQueryLeft is the number of query bases that will extend
    # to the left of the start of the subject in our plots.
    unmatchedQueryLeft = queryStart

    # Set the query offsets into the subject.
    queryStartInSubject = subjectStart - unmatchedQueryLeft
    queryEndInSubject = queryStartInSubject + queryLen + queryGaps

    _sanityCheck(
        subjectStart,
        subjectEnd,
        queryStart,
        queryEnd,
        queryStartInSubject,
        queryEndInSubject,
        hsp,
        queryLen,
        subjectGaps,
        queryGaps,
        locals(),
    )

    return {
        "readStart": queryStart,
        "readEnd": queryEnd,
        "readStartInSubject": queryStartInSubject,
        "readEndInSubject": queryEndInSubject,
        "subjectStart": subjectStart,
        "subjectEnd": subjectEnd,
    }

def printHSP(hsp, indent=''):
    for attr in ['bits', 'expect', 'frame', 'query_end', 'query_start',
                 'sbjct', 'query', 'sbjct_end', 'sbjct_start']:
        print('%s%s: %s' % (indent, attr, hsp[attr]))


def normalizeHSP(hsp, readLen, blastApplication):
    """
    Examine an HSP and return information about where the query and subject
    match begins and ends.  Return a dict with keys that allow the query to
    be displayed against the subject. The returned readStartInSubject and
    readEndInSubject indices are offsets into the subject. I.e., they
    indicate where in the subject the query falls.

    In the returned object, all indices are suitable for Python string
    slicing etc.  We must be careful to convert from the 1-based offsets
    found in BLAST output properly.

    hsp['frame'] is a (query, subject) 2-tuple, with both values coming from
    {-3, -2, -1, 1, 2, 3}. The sign indicates negative or positive sense
    (i.e., the direction of reading through the query or subject to get the
    alignment). The value is the nucleotide match offset modulo 3, plus one
    (i.e., it tells us which of the 3 possible reading frames is used in
    the match). The value is redundant because that information could also
    be obtained from the mod 3 value of the match offset.

    NOTE: the returned readStartInSubject value may be negative.  We consider
    the hit sequence to start at offset 0.  So if the read string has
    sufficient additional nucleotides before the start of the alignment
    match, it may protrude to the left of the hit. Similarly, the returned
    readEndInSubject can be greater than the subjectEnd.

    @param hsp: an HSP in the form of a C{dict}, built from a BLAST record.
        All passed hsp offsets are 1-based.
    @param readLen: the length of the read sequence.
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
        print('normalizeHSP error:')
        print('  readLen: %d' % readLen)
        for var in sorted(locals):
            if var in ('debugPrint', 'hsp'):
                continue
            print('  %s: %s' % (var, locals[var]))
        print('  Original HSP:')
        printHSP(hsp, '    ')
        if msg:
            raise AssertionError(msg)
        else:
            raise AssertionError()

    readPositive = hsp['frame'][0] > 0
    hitPositive = hsp['frame'][1] > 0

    # The following variable names with underscores match the names of
    # attributes BioPython uses and the values (1-based) match those
    # reported by BLAST.
    read_start = hsp['query_start']
    read_end = hsp['query_end']
    sbjct_start = hsp['sbjct_start']
    sbjct_end = hsp['sbjct_end']

    # When the read is positive, BLASTN and TBLASTX give read offsets
    # ascending.
    #
    # TBLASTX reports negative read sense with indices ascending.
    # BLASTN does not report negative read sense.
    #
    # In all cases the read offsets should be ascending.
    if read_start > read_end:
        debugPrint(locals(),
                   'Assertion "read_start <= read_end" failed. Read '
                   'positive is %s. read_start = %d, read_end = %d' %
                   (readPositive, read_start, read_end))

    if hitPositive:
        # Make sure indices are ascending.
        if sbjct_start > sbjct_end:
            debugPrint(locals())
    else:
        # Hit is negative. Its indices will be ascending for TBLASTX
        # output but descending for BLASTN :-( Make sure we have them
        # ascending.
        if sbjct_start > sbjct_end:
            sbjct_start, sbjct_end = sbjct_end, sbjct_start

    # Now that we have asserted what we can about the original HSP values
    # and gotten them into ascending order, make some sane 0-based offsets.
    readStartInSubject = read_start - 1
    readEndInSubject = read_end
    subjectStart = sbjct_start - 1
    subjectEnd = sbjct_end

    if blastApplication == 'blastx':
        # In Blastx output, hit offsets are based on protein sequence
        # length but queries (and the reported offsets) are nucleotide.
        # Convert the read offsets to protein because we will plot against
        # the hit (protein).
        #
        # Note that readStartInSubject and readEndInSubject may not be 0 mod
        # 3. They are offsets into the read string giving the position of
        # the AA, which depends on the translation frame.
        readStartInSubject = int(readStartInSubject / 3)
        readEndInSubject = int(readEndInSubject / 3)

    # No operations on original 1-based HSP variables (with underscores)
    # should appear beyond this point.

    subjectLength = subjectEnd - subjectStart
    readLength = readEndInSubject - readStartInSubject

    # NOTE: readLength (above) is a really bad name. It's actually going to
    # hold the length of the match in the query. I don't know why
    # readEndInSubject - readStartInSubject is used (I mean why those two
    # variables are not named readEnd and readStart). Maybe someone made a
    # find and replace editing error which changed their names. Anyway, the
    # readLength variable is confusingly named because this function is
    # passed a 'readLen' argument, which does happen to be the full length
    # of the read.  This should be cleaned up. See ../diamond/hsp.py for
    # something cleaner.

    hitGaps = hsp['sbjct'].count('-')
    readGaps = hsp['query'].count('-')

    # Sanity check that the length of the matches in the hit and read
    # are identical, taking into account gaps in either (indicated by '-'
    # characters in the match sequences, as returned by BLAST).
    subjectLengthWithGaps = subjectLength + hitGaps
    readLengthWithGaps = readLength + readGaps
    if subjectLengthWithGaps != readLengthWithGaps:
        debugPrint(locals(),
                   'Including gaps, hit match length (%d) != Read match '
                   'length (%d)' % (subjectLengthWithGaps,
                                    readLengthWithGaps))

    # TODO: check the mod 3 value of the start offsets.

    # Calculate read indices. These are indices relative to the hit!

    # unmatchedReadLeft is the number of read bases that will be sticking
    # out to the left of the start of the hit in our plots.
    if readPositive:
        unmatchedReadLeft = readStartInSubject
    else:
        unmatchedReadLeft = readLen - readEndInSubject

    # Set the read offsets based on the direction the match with the
    # hit takes.
    if hitPositive:
        readStartInSubject = subjectStart - unmatchedReadLeft
        readEndInSubject = readStartInSubject + readLen + readGaps
    else:
        readEndInSubject = subjectEnd + unmatchedReadLeft

        readStartInSubject = readEndInSubject - readLen - readGaps

    # Final sanity checks.
    if readStartInSubject > subjectStart:
        debugPrint(locals(), 'readStartInSubject > subjectStart')
    if readEndInSubject < subjectEnd:
        debugPrint(locals(), 'readEndInSubject < subjectEnd')

    return {
        'readStart': read_start - 1,
        'readEnd': read_end,
        'readStartInSubject': readStartInSubject,
        'readEndInSubject': readEndInSubject,
        'subjectStart': subjectStart,
        'subjectEnd': subjectEnd,
    }

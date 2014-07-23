from math import exp
from functools import total_ordering

from dark.hsp import HSP


@total_ordering
class EValueHSP(HSP):
    """
    Holds information about a high-scoring pair from a read alignment.
    Comparisons are done as for e-values (smaller is better).

    All offsets are zero-based and follow the Python convention that
    the 'end' positions are not included in the string.

    See dark.hsp for the __init__ implementation.

    @param readStart: The offset in the read where the match begins.
    @param readEnd: The offset in the read where the match ends.
    @param readStartInHit: The offset in the hit where the match of
        the read starts.
    @param readEndInHit: The offset in the hit where the match of
        the read ends.
    @param hitStart: The offset in the hit where the match begins.
    @param hitEnd: The offset in the hit where the match ends.
    @param readMatchedSequence: The matched part of the read. Note that
        this may contain gaps (marked with '-').
    @param hitMatchedSequence: The matched part of the hit. Note that
        this may contain gaps (marked with '-').
    @param score: The score for the hit. Large values are considered
        better.
    """

    def __lt__(self, other):
        return self.score > other.score

    def __eq__(self, other):
        return self.score == other.score

    def betterThan(self, score):
        """
        Compare this HSP's score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this HSP's score is the better.
        """
        return self.score < score


def printHSP(hsp, indent=''):
    for attr in ['bits', 'expect', 'frame', 'query_end', 'query_start',
                 'sbjct', 'query', 'sbjct_end', 'sbjct_start']:
        print '%s%s: %s' % (indent, attr, hsp[attr])


def normalizeHSP(hsp, readLen, blastApplication):
    """
    Examine the sense of an HSP and return information about where the
    read and the alignment (match) begin and end.  Return a dict with keys
    that allow the read and the alignment to be displayed relative to the
    hit orientation (i.e., with start < stop for both the read and the
    match). The returned read indices are offsets into the hit. I.e.,
    they indicate where on the hit the read lies.

    In the returned object, all indices are suitable for Python string
    slicing etc.  We must be careful to convert from the 1-based offsets
    found in BLAST output properly.

    hsp['frame'] is a (read, hit) 2-tuple, with both values coming from
    {-3, -2, -1, 1, 2, 3}. The sign indicates negative or positive sense
    (i.e., the direction of reading through the read or hit to get the
    alignment). The value is the nucleotide match offset modulo 3, plus one
    (i.e., it tells us which of the 3 possible reading frames is used in
    the match). The value is redundant because that information could also
    be obtained from the mod 3 value of the match offset.

    NOTE: the returned readStartInHit value may be negative.  We consider
    the hit sequence to start at offset 0.  So if the read string has
    sufficient additional nucleotides before the start of the alignment
    match, it may protrude to the left of the hit. Similarly, the returned
    readEndInHit can be greater than the hitEnd.

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
        print 'normalizeHSP error:'
        print '  readLen: %d' % readLen
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
    readStartInHit = read_start - 1
    readEndInHit = read_end
    hitStart = sbjct_start - 1
    hitEnd = sbjct_end

    if blastApplication == 'blastx':
        # In Blastx output, hit offsets are based on protein sequence
        # length but queries (and the reported offsets) are nucleotide.
        # Convert the read offsets to protein because we will plot against
        # the hit (protein).
        #
        # Note that readStartInHit and readEndInHit may not be 0 mod
        # 3. They are offsets into the read string giving the position of
        # the AA, which depends on the translation frame.
        readStartInHit = int(readStartInHit / 3)
        readEndInHit = int(readEndInHit / 3)

    # No operations on original 1-based HSP variables (with underscores)
    # should appear beyond this point.

    hitLength = hitEnd - hitStart
    readLength = readEndInHit - readStartInHit

    hitGaps = hsp['sbjct'].count('-')
    readGaps = hsp['query'].count('-')

    # Sanity check that the length of the matches in the hit and read
    # are identical, taking into account gaps in either (indicated by '-'
    # characters in the match sequences, as returned by BLAST).
    hitLengthWithGaps = hitLength + hitGaps
    readLengthWithGaps = readLength + readGaps
    if hitLengthWithGaps != readLengthWithGaps:
        debugPrint(locals(),
                   'Including gaps, hit match length (%d) != Read match '
                   'length (%d)' % (hitLengthWithGaps,
                                    readLengthWithGaps))

    # TODO: check the mod 3 value of the start offsets.

    # Calculate read indices. These are indices relative to the hit!

    # unmatchedReadLeft is the number of read bases that will be sticking
    # out to the left of the start of the hit in our plots.
    if readPositive:
        unmatchedReadLeft = readStartInHit
    else:
        unmatchedReadLeft = readLen - readEndInHit

    # Set the read offsets based on the direction the match with the
    # hit takes.
    if hitPositive:
        readStartInHit = hitStart - unmatchedReadLeft
        readEndInHit = readStartInHit + readLen + readGaps
    else:
        readEndInHit = hitEnd + unmatchedReadLeft

        readStartInHit = readEndInHit - readLen - readGaps

    # Final sanity checks.
    if readStartInHit > hitStart:
        debugPrint(locals(), 'readStartInHit > hitStart')
    if readEndInHit < hitEnd:
        debugPrint(locals(), 'readEndInHit < hitEnd')

    return {
        'readStart': read_start - 1,
        'readEnd': read_end,
        'readStartInHit': readStartInHit,
        'readEndInHit': readEndInHit,
        'hitStart': hitStart,
        'hitEnd': hitEnd,
    }

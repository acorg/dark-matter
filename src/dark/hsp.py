from functools import total_ordering

from dark.score import HigherIsBetterScore, LowerIsBetterScore


@total_ordering
class _Base:
    """
    Holds information about a matching region from a read alignment.

    You should not use this class directly. Use one of its subclasses,
    either HSP or LSP, depending on whether you want numerically higher
    scores to be considered better (HSP) or worse (LSP).

    Below is an example alignment to show the locations of the six
    start/end offsets.  The upper four are offsets into the subject. The
    lower two are offsets into the read. Note that the read has two gaps
    ('-' characters). All offsets are zero-based and follow the Python
    convention that the 'end' positions are not included in the string.

                           readStartInSubject              readEndInSubject
                           |                               |
                           |                               |
                           |   subjectStart    subjectEnd  |
                           |   |               |           |
                           |   |               |           |
    Subject:  .................ACGTAAAGGCTTAGGT.................
    Read:                  ....ACGTA-AGGCTT-GGT............
                               |               |
                               |               |
                               readStart       readEnd

    Note that the above is just one alignment, and that others are possible
    (e.g., with the read extending beyond the end(s) of the subject, or the
    subject also with gaps in it). The point of the example diagram is to show
    what the six variable names will always refer to, not to enumerate all
    possible alignments (the tests in test/blast/test_hsp.py go through
    many different cases). The classes in this file are just designed to hold
    the variables associated with an HSP and to make it easy to compare them.

    @param readStart: The offset in the read where the match begins.
    @param readEnd: The offset in the read where the match ends.
    @param readStartInSubject: The offset in the subject where the match of
        the read starts.
    @param readEndInSubject: The offset in the subject where the match of
        the read ends.
    @param readFrame: The reading frame for the read, a value from
        {-3, -2, -1, 1, 2, 3} where the sign indicates negative or positive
        sense.
    @param subjectStart: The offset in the subject where the match begins.
    @param subjectEnd: The offset in the subject where the match ends.
    @param subjectFrame: The reading frame for the subject, a value from
        {-3, -2, -1, 1, 2, 3} where the sign indicates negative or positive
        sense.
    @param readMatchedSequence: The matched part of the read. Note that
        this may contain gaps (marked with '-').
    @param subjectMatchedSequence: The matched part of the subject. Note that
        this may contain gaps (marked with '-').
    @param identicalCount: The C{int} number of positions at which the subject
        and query were identical.
    @param positiveCount: The C{int} number of positions at which the subject
        and query had a positive score in the scoring matrix used during
        matching (this is probably only different from the C{identicalCount}
        when matching amino acids (i.e., not nucleotides).
    @param percentIdentical: A C{float} percentage (i.e., ranging from 0.0 to
        100.0, NOT a fraction) of amino acids that were identical in the match.
    """

    def __init__(
        self,
        readStart=None,
        readEnd=None,
        readStartInSubject=None,
        readEndInSubject=None,
        readFrame=None,
        subjectStart=None,
        subjectEnd=None,
        subjectFrame=None,
        readMatchedSequence=None,
        subjectMatchedSequence=None,
        identicalCount=None,
        percentIdentical=None,
        positiveCount=None,
        percentPositive=None,
    ):
        self.readStart = readStart
        self.readEnd = readEnd
        self.readStartInSubject = readStartInSubject
        self.readEndInSubject = readEndInSubject
        self.readFrame = readFrame
        self.subjectStart = subjectStart
        self.subjectEnd = subjectEnd
        self.subjectFrame = subjectFrame
        self.readMatchedSequence = readMatchedSequence
        self.subjectMatchedSequence = subjectMatchedSequence
        self.identicalCount = identicalCount
        self.percentIdentical = percentIdentical
        self.positiveCount = positiveCount
        self.percentPositive = percentPositive

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.score == other.score

    def betterThan(self, score):
        """
        Compare this instance's score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this score is the better.
        """
        return self.score.betterThan(score)

    def toDict(self):
        """
        Get information about the HSP/LSP as a dictionary.

        @return: A C{dict} representation of the HSP/LSP.
        """
        return {
            "readStart": self.readStart,
            "readEnd": self.readEnd,
            "readStartInSubject": self.readStartInSubject,
            "readEndInSubject": self.readEndInSubject,
            "readFrame": self.readFrame,
            "subjectStart": self.subjectStart,
            "subjectEnd": self.subjectEnd,
            "subjectFrame": self.subjectFrame,
            "readMatchedSequence": self.readMatchedSequence,
            "subjectMatchedSequence": self.subjectMatchedSequence,
            "identicalCount": self.identicalCount,
            "percentIdentical": self.percentIdentical,
            "positiveCount": self.positiveCount,
            "percentPositive": self.percentPositive,
        }


class HSP(_Base):
    """
    Holds information about a high-scoring pair from a read alignment.
    Comparisons are done as for BLAST or DIAMOND bit scores (higher is better).

    @param score: The numeric score of this HSP.
    """

    def __init__(self, score, **kwargs):
        _Base.__init__(self, **kwargs)
        self.score = HigherIsBetterScore(score)

    def toDict(self):
        """
        Get information about the HSP as a dictionary.

        @return: A C{dict} representation of the HSP.
        """
        result = _Base.toDict(self)
        result["score"] = self.score.score
        return result


class LSP(_Base):
    """
    Holds information about a low-scoring pair from a read alignment.
    Comparisons are done as for BLAST or DIAMOND e-values (smaller is better).

    @param score: The numeric score of this LSP.
    """

    def __init__(self, score, **kwargs):
        _Base.__init__(self, **kwargs)
        self.score = LowerIsBetterScore(score)

    def toDict(self):
        """
        Get information about the LSP as a dictionary.

        @return: A C{dict} representation of the LSP.
        """
        result = _Base.toDict(self)
        result["score"] = self.score.score
        return result

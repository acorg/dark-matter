from functools import total_ordering

from dark.score import HigherIsBetterScore, LowerIsBetterScore


@total_ordering
class _Base(object):
    """
    Holds information about a matching region from a read alignment.

    You should not use this class directly. Use one of its subclasses,
    either HSP or LSP, depending on whether you want numerically higher
    scores to be considered better (HSP) or worse (LSP).

    All offsets are zero-based and follow the Python convention that
    the 'end' positions are not included in the string.

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
    """
    def __init__(self, readStart=None, readEnd=None, readStartInHit=None,
                 readEndInHit=None, hitStart=None, hitEnd=None,
                 readMatchedSequence=None, hitMatchedSequence=None):
        self.readStart = readStart
        self.readEnd = readEnd
        self.readStartInHit = readStartInHit
        self.readEndInHit = readEndInHit
        self.hitStart = hitStart
        self.hitEnd = hitEnd
        self.readMatchedSequence = readMatchedSequence
        self.hitMatchedSequence = hitMatchedSequence

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.score == other.score

    def betterThan(self, score):
        """
        Compare this LSP's score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this score is the better.
        """
        return self.score.betterThan(score)


class HSP(_Base):
    """
    Holds information about a high-scoring pair from a read alignment.
    Comparisons are done as for BLAST bit scores (higher is better).

    @param score: The numeric score of this HSP.
    """

    def __init__(self, score, **kwargs):
        _Base.__init__(self, **kwargs)
        self.score = HigherIsBetterScore(score)


class LSP(_Base):
    """
    Holds information about a low-scoring pair from a read alignment.
    Comparisons are done as for BLAST e-values (smaller is better).

    @param score: The numeric score of this LSP.
    """

    def __init__(self, score, **kwargs):
        _Base.__init__(self, **kwargs)
        self.score = LowerIsBetterScore(score)

from functools import total_ordering


@total_ordering
class HSP(object):
    """
    Holds information about a high-scoring pair from a read alignment.

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
    @param score: The score for the hit. Large values are considered
        better.
    """
    def __init__(self, readStart=None, readEnd=None,
                 readStartInHit=None, readEndInHit=None,
                 hitStart=None, hitEnd=None,
                 readMatchedSequence=None, hitMatchedSequence=None,
                 score=None):
        self.readStart = readStart
        self.readEnd = readEnd
        self.readStartInHit = readStartInHit
        self.readEndInHit = readEndInHit
        self.hitStart = hitStart
        self.hitEnd = hitEnd
        self.readMatchedSequence = readMatchedSequence
        self.hitMatchedSequence = hitMatchedSequence
        self.score = score

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.score == other.score

    def betterThan(self, score):
        """
        Compare this HSP's score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this HSP's score is the better.
        """
        return self.score > score

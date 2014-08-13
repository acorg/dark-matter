from functools import total_ordering


@total_ordering
class HigherIsBetterScore(object):
    """
    Provide comparison functions for scores where numerically higher values
    are considered better.

    @param score: The numeric score of this HSP.
    """
    def __init__(self, score):
        self.score = score

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.score == other.score

    def betterThan(self, score):
        """
        Compare this score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this score is the better.
        """
        return self.score > score


@total_ordering
class LowerIsBetterScore(object):
    """
    Provide comparison functions for scores where numerically lower values
    are considered better.

    @param score: The numeric score of this LSP.
    """
    def __init__(self, score):
        self.score = score

    def __lt__(self, other):
        return self.score > other.score

    def __eq__(self, other):
        return self.score == other.score

    def betterThan(self, score):
        """
        Compare this score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this score is the better.
        """
        return self.score < score

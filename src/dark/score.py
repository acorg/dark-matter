from functools import total_ordering


@total_ordering
class HigherIsBetterScore:
    """
    Provide comparison functions for scores where numerically higher values
    are considered better.

    @param score: The numeric score of this HSP.
    """

    def __init__(self, score: float) -> None:
        self.score = score

    def __lt__(self, other: object) -> bool:
        if isinstance(other, HigherIsBetterScore):
            return self.score < other.score
        return NotImplemented

    def __eq__(self, other: object) -> bool:
        if isinstance(other, HigherIsBetterScore):
            return self.score == other.score
        return NotImplemented

    def betterThan(self, score: float) -> bool:
        """
        Compare this score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this score is the better.
        """
        return self.score > score


@total_ordering
class LowerIsBetterScore:
    """
    Provide comparison functions for scores where numerically lower values
    are considered better.

    @param score: The numeric score of this LSP.
    """

    def __init__(self, score: float) -> None:
        self.score = score

    def __lt__(self, other: object) -> bool:
        if isinstance(other, LowerIsBetterScore):
            return self.score > other.score
        return NotImplemented

    def __eq__(self, other: object) -> bool:
        if isinstance(other, LowerIsBetterScore):
            return self.score == other.score
        return NotImplemented

    def betterThan(self, score: float) -> bool:
        """
        Compare this score with another score.

        @param score: A C{float} score.
        @return: A C{bool}, C{True} if this score is the better.
        """
        return self.score < score

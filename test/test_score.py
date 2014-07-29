from unittest import TestCase

from dark.score import HigherIsBetterScore, LowerIsBetterScore


class TestHigherIsBetterScore(TestCase):
    """
    Tests of the C{dark.score.HigherIsBetterScore} class.
    """

    def testExpectedAttr(self):
        """
        The expected score attribute must be present.
        """
        score = HigherIsBetterScore(10)
        self.assertEqual(10, score.score)

    def testEqual(self):
        """
        __eq__ must work as expected.
        """
        self.assertEqual(HigherIsBetterScore(10), HigherIsBetterScore(10))

    def testLessThan(self):
        """
        __lt__ must work as expected.
        """
        self.assertTrue(HigherIsBetterScore(10) < HigherIsBetterScore(20))

    def testBetterThan(self):
        """
        betterThan must work as expected.
        """
        self.assertTrue(HigherIsBetterScore(10).betterThan(5))


class TestLowerIsBetterScore(TestCase):
    """
    Tests of the C{dark.score.LowerIsBetterScore} class.
    """

    def testExpectedAttr(self):
        """
        The expected score attribute must be present.
        """
        score = LowerIsBetterScore(10)
        self.assertEqual(10, score.score)

    def testEqual(self):
        """
        __eq__ must work as expected.
        """
        self.assertEqual(LowerIsBetterScore(10), LowerIsBetterScore(10))

    def testLessThan(self):
        """
        __lt__ must work as expected.
        """
        self.assertTrue(LowerIsBetterScore(20) < LowerIsBetterScore(10))

    def testBetterThan(self):
        """
        betterThan must work as expected.
        """
        self.assertTrue(LowerIsBetterScore(10).betterThan(50))

from unittest import TestCase

from dark.colors import ColorsForCounts


class TestColorsForCounts(TestCase):
    """
    Test the ColorsForCounts class.
    """

    def testCountNotInt(self):
        """
        A C{ValueError} must be raised if a non-integer count is given.
        """
        error = (
            r"^color arguments must be given as space-separated pairs of "
            r"\"count color\" where the count is an integer threshold\. "
            r"Your value \('non-int'\) was not an integer\.$"
        )
        self.assertRaisesRegex(ValueError, error, ColorsForCounts, ("non-int red",))

    def testCountNegative(self):
        """
        A C{ValueError} must be raised if a count is less than zero.
        """
        error = (
            r"^color arguments must be given as space-separated pairs "
            r"of \"count color\" where the count is non-negative. Your "
            r"value \(-4\) is less than 0\.$"
        )
        self.assertRaisesRegex(ValueError, error, ColorsForCounts, ("-4 red",))

    def testNoSpace(self):
        """
        A C{ValueError} must be raised if a count/color value has no space.
        """
        error = (
            r"^color arguments must be given as space-separated pairs "
            r"of \"value color\"\. Your value \('4red'\) does not "
            r"contain a space\.$"
        )
        self.assertRaisesRegex(ValueError, error, ColorsForCounts, ("4red",))

    def testRepeatedCount(self):
        """
        A C{ValueError} must be raised if a count is repeated.
        """
        error = r"^repeated color argument count \(4\)\.$"
        self.assertRaisesRegex(ValueError, error, ColorsForCounts, ("4 red", "4 black"))

    def testThresholdForCountNegative(self):
        """
        The thresholdForCount method must raise AssertionError if the count
        is negative.
        """
        colors = ColorsForCounts(("4 red",))
        error = r"^Count \(-4\) cannot be negative\.$"
        self.assertRaisesRegex(AssertionError, error, colors.thresholdForCount, -4)

    def testThresholdForTooLowCount(self):
        """
        The thresholdForCount method must work correctly when the lowest
        threshold is not reached.
        """
        colors = ColorsForCounts(("4 red",))
        self.assertEqual(0, colors.thresholdForCount(2))

    def testThresholdForCountOneThreshold(self):
        """
        The thresholdForCount method must work correctly when there is one
        threshold.
        """
        colors = ColorsForCounts(("4 red",))
        self.assertEqual(4, colors.thresholdForCount(4))
        self.assertEqual(4, colors.thresholdForCount(10))

    def testThresholdForCountTwoThresholds(self):
        """
        The thresholdForCount method must work correctly when there are two
        thresholds.
        """
        colors = ColorsForCounts(("4 red", "10 green"))
        self.assertEqual(4, colors.thresholdForCount(6))
        self.assertEqual(10, colors.thresholdForCount(12))

    def testThresholdForCountThreeThresholds(self):
        """
        The thresholdForCount method must work correctly when there are three
        thresholds.
        """
        colors = ColorsForCounts(("4 red", "10 green", "100 blue"))
        self.assertEqual(4, colors.thresholdForCount(6))
        self.assertEqual(10, colors.thresholdForCount(86))
        self.assertEqual(100, colors.thresholdForCount(1000))

    def testColorForTooLowCount(self):
        """
        The colorForCount method must work correctly when the lowest
        threshold is not reached.
        """
        colors = ColorsForCounts(("4 red",))
        self.assertEqual("black", colors.colorForCount(2))

    def testColorForTooLowCountDefault(self):
        """
        The colorForCount method must work correctly when the lowest
        threshold is not reached and a default color is given.
        """
        colors = ColorsForCounts(("4 red",), defaultColor="fuschia")
        self.assertEqual("fuschia", colors.colorForCount(2))

    def testColorForCountOneThreshold(self):
        """
        The colorForCount method must work correctly when there is one
        threshold.
        """
        colors = ColorsForCounts(("4 red",))
        self.assertEqual("black", colors.colorForCount(1))
        self.assertEqual("red", colors.colorForCount(4))
        self.assertEqual("red", colors.colorForCount(10))

    def testColorForCountTwoThresholds(self):
        """
        The colorForCount method must work correctly when there are two
        thresholds.
        """
        colors = ColorsForCounts(("4 red", "10 green"))
        self.assertEqual("black", colors.colorForCount(1))
        self.assertEqual("red", colors.colorForCount(6))
        self.assertEqual("green", colors.colorForCount(12))

    def testColorForCountThreeThresholds(self):
        """
        The colorForCount method must work correctly when there are three
        thresholds.
        """
        colors = ColorsForCounts(("4 red", "10 green", "100 blue"))
        self.assertEqual("black", colors.colorForCount(1))
        self.assertEqual("red", colors.colorForCount(6))
        self.assertEqual("green", colors.colorForCount(86))
        self.assertEqual("blue", colors.colorForCount(1000))

    def testThresholdToCssName(self):
        """
        The thresholdToCssName method must work as expected.
        """
        colors = ColorsForCounts(("4 red",))
        self.assertEqual("threshold-4", colors.thresholdToCssName(4))

    def testColorsForThreeThresholds(self):
        """
        The colors attribute must be as expected when there are three
        thresholds.
        """
        colors = ColorsForCounts(("4 red", "10 green", "100 blue"), defaultColor="pink")
        self.assertEqual(
            ((100, "blue"), (10, "green"), (4, "red"), (0, "pink")), colors.colors
        )

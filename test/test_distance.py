from unittest import TestCase

from dark.distance import levenshtein


class TestLevenshtein(TestCase):
    """
    Tests for the dark.simplify.levenshtein function.
    """

    def testIdentical(self):
        """
        Two identical strings must have distance zero.
        """
        self.assertEqual(0, levenshtein("BLAH", "BLAH"))

    def testMutation(self):
        """
        Test a single character results in a distance of 1.
        """
        self.assertEqual(1, levenshtein("ACGTACACACG", "ACGTACACACT"))

    def testInsert(self):
        """
        Test a string insertion that results in a distance of 2.
        """
        self.assertEqual(2, levenshtein("AGTACACACTG", "ACGTACACACT"))

from six import assertRaisesRegex
from unittest import TestCase

from dark.btop import countGaps, parseBtop


class TestParseBtop(TestCase):
    """
    Tests for the parseBtop function.
    """
    def testOneLetter(self):
        """
        An argument with just one letter must result in a ValueError.
        """
        error = ("^btop string 'F' has a trailing query letter 'F' with no "
                 "corresponding subject letter$")
        assertRaisesRegex(self, ValueError, error, list, parseBtop('F'))

    def testConsecutiveGaps(self):
        """
        An argument that has two consecutive gaps characters must result in a
        ValueError.
        """
        error = "^btop string '36--' has two consecutive gaps at offset 2$"
        assertRaisesRegex(self, ValueError, error, list, parseBtop('36--'))

    def testConsecutiveIdentical(self):
        """
        An argument that has two consecutive identical (non-gap) characters
        must result in a ValueError.
        """
        error = ("^btop string '36AA' has two consecutive identical 'A' "
                 "letters at offset 2$")
        assertRaisesRegex(self, ValueError, error, list, parseBtop('36AA'))

    def testEmpty(self):
        """
        An empty argument must result in an empty list.
        """
        self.assertEqual([], list(parseBtop('')))

    def testOneNumberWithTrailingOneLetter(self):
        """
        An argument that is a number with a single letter must result in a
        ValueError.
        """
        error = ("^btop string '36F' has a trailing query letter 'F' with no "
                 "corresponding subject letter$")
        assertRaisesRegex(self, ValueError, error, list, parseBtop('36F'))

    def testThreeLetters(self):
        """
        An argument that has three letters must result in a ValueError.
        """
        error = ("^btop string 'ABC' has a trailing query letter 'C' with no "
                 "corresponding subject letter$")
        assertRaisesRegex(self, ValueError, error, list, parseBtop('ABC'))

    def testOneLetterThenANumber(self):
        """
        An argument that is a single letter followed by a number must result
        in a ValueError.
        """
        error = ("^btop string 'F36' has a query letter 'F' at offset 0 with "
                 "no corresponding subject letter$")
        assertRaisesRegex(self, ValueError, error, list, parseBtop('F36'))

    def testTwoNumbersWithOneLetterBetween(self):
        """
        An argument that is a number, a single letter, and another number must
        result in a ValueError.
        """
        error = ("^btop string '36F77' has a query letter 'F' at offset 2 "
                 "with no corresponding subject letter$")
        assertRaisesRegex(self, ValueError, error, list, parseBtop('36F77'))

    def testOneNumber(self):
        """
        An argument that is just one number must give the expected result.
        """
        self.assertEqual([54], list(parseBtop('54')))

    def testOneNumberThatIsZero(self):
        """
        An argument that is just the number zero must give the expected result.
        """
        self.assertEqual([0], list(parseBtop('0')))

    def testOneNumberWithLeadingZeroes(self):
        """
        An argument that is just one number with leading zeroes must give the
        expected result.
        """
        self.assertEqual([54], list(parseBtop('0054')))

    def testOneQuerySubjectPair(self):
        """
        An argument that is a single query/subject letter pair must give the
        expected result.
        """
        self.assertEqual([('A', 'G')], list(parseBtop('AG')))

    def testTwoQuerySubjectPairs(self):
        """
        An argument that has two query/subject letter pairs must give the
        expected result.
        """
        self.assertEqual([('A', 'G'), ('C', 'T')], list(parseBtop('AGCT')))

    def testOneQuerySubjectPairAndANumber(self):
        """
        An argument that is a single query/subject letter pair followed by a
        number must give the expected result.
        """
        self.assertEqual([('A', 'G'), 33], list(parseBtop('AG33')))


class TestCountGaps(TestCase):
    """
    Tests for the countGaps function.
    """
    def testEmpty(self):
        """
        An argument with an empty string must produce the expected result.
        """
        self.assertEqual((0, 0), countGaps(''))

    def testNumberOnly(self):
        """
        An argument with just a number must produce the expected result.
        """
        self.assertEqual((0, 0), countGaps('88'))

    def testLettersButNoGaps(self):
        """
        An argument with just letters must produce the expected result.
        """
        self.assertEqual((0, 0), countGaps('FGAC'))

    def testOneQueryGap(self):
        """
        An argument with just a query gap must produce the expected result.
        """
        self.assertEqual((1, 0), countGaps('-G'))

    def testOneSubjectGap(self):
        """
        An argument with just a subject gap must produce the expected result.
        """
        self.assertEqual((0, 1), countGaps('G-'))

    def testOneQueryAndOneSubjectGap(self):
        """
        An argument with a query and a subject gap must produce the expected
        result.
        """
        self.assertEqual((1, 1), countGaps('G--G'))

    def testMultipleQueryAndSubjectGaps(self):
        """
        An argument with multiple query and a subject gaps must produce the
        expected result.
        """
        self.assertEqual((3, 2), countGaps('-GG-34-T-T39F-'))

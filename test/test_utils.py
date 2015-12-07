import six
from unittest import TestCase

from dark.utils import numericallySortFilenames, median


class TestNumericallySortFilenames(TestCase):
    """
    Test the numericallySortFilenames function.
    """

    def testNoNames(self):
        """
        An empty list must be returned when an empty list is given.
        """
        self.assertEqual([], numericallySortFilenames([]))

    def testOneNonNumericName(self):
        """
        A list with a single non-numeric name should result in that same
        name being returned.
        """
        self.assertEqual(['hey'], numericallySortFilenames(['hey']))

    def testOneNumericName(self):
        """
        A list with a single numeric name should result in that same
        name being returned.
        """
        self.assertEqual(['3.json'], numericallySortFilenames(['3.json']))

    def testSeveralNames(self):
        """
        A list with several numeric names should result in a correctly
        sorted list of names being returned.
        """
        self.assertEqual(
            ['1.json', '2.json', '3.json'],
            numericallySortFilenames(['3.json', '1.json', '2.json']))

    def testSeveralNamesWithUnequalPrefixLengths(self):
        """
        A list with several numeric names whose numeric prefixes differ
        in length should result in a correctly sorted list of names being
        returned.
        """
        self.assertEqual(
            ['2.json', '3.json', '21.json', '35.json', '250.json'],
            numericallySortFilenames(
                ['3.json', '21.json', '35.json', '250.json', '2.json']))

    def testBasename(self):
        """
        Sorting must be according to file basename.
        """
        self.assertEqual(
            ['../output/2.json', '../output/3.json', '../output/21.json',
             '../output/35.json', '../output/250.json'],
            numericallySortFilenames(
                ['../output/3.json', '../output/21.json', '../output/35.json',
                 '../output/250.json', '../output/2.json']))


class TestMedian(TestCase):
    """
    Tests for the median function.
    """

    def testEmptyArgRaises(self):
        """
        An empty list must cause median to raise ValueError.
        """
        error = '^arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error, median, [])

    def testMedianOfOne(self):
        """
        The median function must work on a list of length one.
        """
        self.assertEqual(3, median([3]))

    def testMedianOfTwo(self):
        """
        The median function must work on a list of length two.
        """
        self.assertEqual(4.5, median([3.1, 5.9]))

    def testMedianOfThree(self):
        """
        The median function must work on a list of length threee.
        """
        self.assertEqual(5.9, median([3.1, 7.6, 5.9]))

    def testMedianOfFour(self):
        """
        The median function must work on a list of length four.
        """
        self.assertEqual(4.5, median([3.1, 1.3, 7.6, 5.9]))

    def testMedianOfFive(self):
        """
        The median function must work on a list of length five.
        """
        self.assertEqual(5.9, median([3.1, 1.3, 7.6, 9.9, 5.9]))

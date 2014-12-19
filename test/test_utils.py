from unittest import TestCase

from dark.utils import numericallySortFilenames


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

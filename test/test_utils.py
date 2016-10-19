import six
import bz2
import gzip
from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import mockOpen

from dark.utils import numericallySortFilenames, median, asHandle


class BZ2(object):
    """
    A BZ2File mock.
    """
    def __init__(self, data):
        self._data = data

    def read(self):
        return self._data

GZIP = BZ2


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


class TestAsHandle(TestCase):
    """
    Test the asHandle function
    """

    def testOpenFile(self):
        """
        When an open file pointer is passed to asHandle, that same file
        pointer must be returned.
        """
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            fp = open('file')
            with asHandle(fp) as newfp:
                self.assertIs(fp, newfp)

    def testStr(self):
        """
        When a string filename is passed to asHandle, it must be possible to
        read the correct data from the fp that is returned.
        """
        mockOpener = mockOpen(read_data='xxx')
        with patch.object(builtins, 'open', mockOpener):
            with asHandle('file') as fp:
                self.assertEqual('xxx', fp.read())

    def testBZ2(self):
        """
        When a string '*.bz2' filename is passed to asHandle, it must be
        possible to read the correct data from the fp that is returned.
        """
        # This test should be better. It should actually create some bz2
        # compressed data and make sure that it's decompressed
        # properly. But Python mocking makes me so confused...
        result = BZ2('xxx')

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.return_value = result
            with asHandle('file.bz2') as fp:
                self.assertEqual('xxx', fp.read())

    def testGzip(self):
        """
        When a string '*.gz' filename is passed to asHandle, it must be
        possible to read the correct data from the fp that is returned.
        """
        # This test should be better. It should actually create some gzip
        # compressed data and make sure that it's decompressed
        # properly. But Python mocking makes me so confused...
        result = GZIP('xxx')

        with patch.object(gzip, 'GzipFile') as mockMethod:
            mockMethod.return_value = result
            with asHandle('file.gz') as fp:
                self.assertEqual('xxx', fp.read())

from six import assertRaisesRegex
from unittest import TestCase
from tempfile import mkstemp
from os import close, unlink, write
from contextlib import contextmanager

from dark.reads import Read
from dark.sam import PaddedSAM, UnequalReferenceLengthError, UnknownReference


# These tests actually use the filesystem to read files. That's due to the API
# to pysam and the fact that it calls a C function to open files, so we can't
# mock Python's 'open' method. Hence the following context manager.
@contextmanager
def dataFile(data):
    """
    Create a context manager to store data in a temporary file and
    later remove it.
    """
    fd, filename = mkstemp()
    write(fd, data.encode('utf-8'))
    close(fd)
    yield filename
    unlink(filename)


class TestPaddedSAM(TestCase):
    """
    Test the PaddedSAM class.
    """
    def testReferencesToStr(self):
        """
        The referencesToStr method must return the expected string.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            self.assertEqual('id1 (length 90)\nid2 (length 91)',
                             ps.referencesToStr())
            ps.close()

    def testUnknownReferences(self):
        """
        Passing an unknown reference name to 'queries' must result in an
        UnknownReference exception.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            error = ("^Reference 'unknown' is not present in the "
                     "SAM/BAM file\\.$")
            queries = ps.queries(referenceName='unknown')
            assertRaisesRegex(self, UnknownReference, error, list, queries)
            ps.close()

    def testUnequalReferenceLengths(self):
        """
        Passing no reference name to 'queries' when the references have
        different lengths must result in an UnequalReferenceLengthError
        exception.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            error = ('^Your SAM/BAM file has 2 reference sequences, and their '
                     'lengths \(90, 91\) are not all identical\.$')
            queries = ps.queries()
            assertRaisesRegex(self, UnequalReferenceLengthError, error, list,
                              queries)
            ps.close()

    def testSimpleMatch(self):
        """
        A simple all-'M' match must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 6M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testRcNeeded(self):
        """
        A reverse-complimented match (flag = 16) when rcNeeded=True is passed
        must result in the expected (reverse complimented) padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 16 ref1 2 60 6M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(rcNeeded=True))
            self.assertEqual(Read('query1', '-CCTAGA---'), read)
            ps.close()

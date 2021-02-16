from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch, mock_open
except ImportError:
    from mock import patch

from dark.filter import ReadSetFilter, TitleFilter
from dark.reads import Read
from dark.titles import TitleAlignment, TitleAlignments


class TitleFilterTest(TestCase):
    """
    Tests for the L{dark.filter.TitleFilter} class.
    """

    def testNoRestriction(self):
        """
        Testing for acceptance against a title filter that has no
        restrictions should return C{TitleFilter.DEFAULT_ACCEPT}.
        """
        tf = TitleFilter()
        self.assertEqual(TitleFilter.DEFAULT_ACCEPT, tf.accept('hey'))

    def testPositiveRegex(self):
        """
        Testing for acceptance against a title filter with a positive regex
        must work.
        """
        tf = TitleFilter(positiveRegex=r'x+\s')
        self.assertEqual(TitleFilter.DEFAULT_ACCEPT, tf.accept('hey xxx you'))
        self.assertEqual(TitleFilter.REJECT, tf.accept('hey xxyou'))

    def testNegativeRegex(self):
        """
        Testing for acceptance against a title filter with a negative regex
        must work.
        """
        tf = TitleFilter(negativeRegex=r'x+\s')
        self.assertEqual(TitleFilter.REJECT, tf.accept('hey xxx you'))
        self.assertEqual(TitleFilter.DEFAULT_ACCEPT, tf.accept('hey xxyou'))

    def testPositiveRegexHasPrecedenceOverRepeatedTruncatedTitle(self):
        """
        Testing for acceptance against a title filter with a positive regex
        must have precedence over checking for truncated titles when the same
        non-matching title (that will be truncated) is passed twice.
        """
        tf = TitleFilter(positiveRegex=r'xxxxx', truncateAfter='virus')
        self.assertEqual(TitleFilter.REJECT, tf.accept('spotty virus 1'))
        self.assertEqual(TitleFilter.REJECT, tf.accept('spotty virus 1'))

    def testNegativeRegexHasPrecedenceOverRepeatedTruncatedTitle(self):
        """
        Testing for acceptance against a title filter with a negative regex
        must have precedence over checking for truncated titles when the same
        matching title (that will be truncated) is passed twice.
        """
        tf = TitleFilter(negativeRegex=r'spotty', truncateAfter='virus')
        self.assertEqual(TitleFilter.REJECT, tf.accept('spotty virus 1'))
        self.assertEqual(TitleFilter.REJECT, tf.accept('spotty virus 1'))

    def testFullWordTruncation(self):
        """
        Testing for acceptance against a title filter with title truncation
        in effect must work if the title contains the C{truncateAfter} string
        as a distint word.
        """
        tf = TitleFilter(truncateAfter=r'virus')
        # Note that the truncation code will chop off the first part of the
        # title (the title ID).
        self.assertEqual(TitleFilter.DEFAULT_ACCEPT,
                         tf.accept('gi|400684|gb|AY421767.1| herpes virus 1'))
        self.assertEqual(TitleFilter.REJECT,
                         tf.accept('gi|400684|gb|AY421767.1| herpes virus 2'))

    def testPartialWordTruncation(self):
        """
        Testing for acceptance against a title filter with title truncation
        in effect must work if the title contains the C{truncateAfter} string
        as a partial word.
        """
        tf = TitleFilter(truncateAfter=r'virus')
        # Note that the truncation code will chop off the first part of the
        # title (the title ID).
        self.assertEqual(TitleFilter.DEFAULT_ACCEPT,
                         tf.accept('gi|400684|gb|AY421767.1| rotavirus 1'))
        self.assertEqual(TitleFilter.REJECT,
                         tf.accept('gi|400684|gb|AY421767.1| rotavirus 2'))

    def testWordTruncationRepeat(self):
        """
        Testing for acceptance against a title filter with title truncation
        in effect must allow the exact same title twice, even if the title
        is being truncated.
        """
        tf = TitleFilter(truncateAfter=r'virus')
        # Note that the truncation code will chop off the first part of the
        # title (the title ID).
        self.assertEqual(TitleFilter.DEFAULT_ACCEPT,
                         tf.accept('gi|400684|gb|AY421767.1| herpes virus 1'))
        self.assertEqual(TitleFilter.DEFAULT_ACCEPT,
                         tf.accept('gi|400684|gb|AY421767.1| herpes virus 1'))

    def testWhitelist(self):
        """
        Testing for acceptance against a title filter with a whitelist
        must work even when a title is ruled out for other violations.
        """
        tf = TitleFilter(whitelist=['always ok'], negativeRegex='ok')
        self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('always ok'))
        self.assertEqual(TitleFilter.REJECT, tf.accept('always ok not'))

    def testBlacklist(self):
        """
        Testing for acceptance against a title filter with a blacklist
        must work.
        """
        tf = TitleFilter(blacklist=['never ok'], positiveRegex='ok')
        self.assertEqual(TitleFilter.REJECT, tf.accept('never ok'))

    def testBlacklistFile(self):
        """
        Testing for acceptance against a title filter with a blacklist file.
        """
        data = '\n'.join(['id1', 'id2']) + '\n'
        with patch.object(builtins, 'open', mock_open(read_data=data)):
            tf = TitleFilter(blacklistFile='black.txt')
            self.assertEqual(TitleFilter.REJECT, tf.accept('id1'))
            self.assertEqual(TitleFilter.REJECT, tf.accept('id2'))
            self.assertEqual(TitleFilter.DEFAULT_ACCEPT, tf.accept('id3'))

    def testBlacklistFileAndBlacklist(self):
        """
        Testing for acceptance against a title filter with a blacklist file and
        some specific other blacklist titles.
        """
        data = '\n'.join(['id1', 'id2']) + '\n'
        with patch.object(builtins, 'open', mock_open(read_data=data)):
            tf = TitleFilter(blacklistFile='black.txt', blacklist=set(['id3']))
            self.assertEqual(TitleFilter.REJECT, tf.accept('id1'))
            self.assertEqual(TitleFilter.REJECT, tf.accept('id2'))
            self.assertEqual(TitleFilter.REJECT, tf.accept('id3'))
            self.assertEqual(TitleFilter.DEFAULT_ACCEPT, tf.accept('id4'))

    def testWhitelistTakesPrecedenceOverBlacklist(self):
        """
        Testing for acceptance against a title filter with a whitelist
        and a blacklist that contain the same title must work (the whitelist
        takes precedence).
        """
        tf = TitleFilter(whitelist=['always ok'], blacklist=['always ok'])
        self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('always ok'))

    def testWhitelistOnly(self):
        """
        Testing for acceptance against a title filter with a whitelist
        and a negative regex that matches everything.
        """
        tf = TitleFilter(whitelist=['always ok'], negativeRegex='.')
        self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('always ok'))
        self.assertEqual(TitleFilter.REJECT, tf.accept('always not ok'))
        self.assertEqual(TitleFilter.REJECT, tf.accept('rubbish'))

    def testWhitelistFileOnly(self):
        """
        Testing for acceptance against a title filter with a whitelist file
        and a negative regex that matches everything.
        """
        data = '\n'.join(['id1', 'id2']) + '\n'
        with patch.object(builtins, 'open', mock_open(read_data=data)):
            tf = TitleFilter(whitelistFile='white.txt', negativeRegex='.')
            self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('id1'))
            self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('id2'))
            self.assertEqual(TitleFilter.REJECT, tf.accept('id3'))

    def testWhitelistFileAndWhitelistOnly(self):
        """
        Testing for acceptance against a title filter with a whitelist file
        and some specific whitelist titles, with a negative regex that matches
        everything.
        """
        data = '\n'.join(['id1', 'id2']) + '\n'
        with patch.object(builtins, 'open', mock_open(read_data=data)):
            tf = TitleFilter(whitelistFile='white.txt', whitelist=set(['id3']),
                             negativeRegex='.')
            self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('id1'))
            self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('id2'))
            self.assertEqual(TitleFilter.WHITELIST_ACCEPT, tf.accept('id3'))
            self.assertEqual(TitleFilter.REJECT, tf.accept('id4'))


class ReadSetTest(TestCase):
    """
    Tests for the L{dark.filter.ReadSetFilter} class.
    """

    def makeTitleAlignments(self, *readIds):
        """
        Create a TitleAlignments instance containing reads with the
        ids given by C{ids}.

        param readIds: A C{list} of integer ids for reads.
        @return: A C{TitleAlignments} instance with reads with the given ids.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        for readId in readIds:
            titleAlignment = TitleAlignment(Read('id' + str(readId), 'A'), [])
            titleAlignments.addAlignment(titleAlignment)
        return titleAlignments

    def testFirstUse(self):
        """
        Testing for acceptance against a read set filter that has not been
        used should return C{True}.
        """
        titleAlignments = self.makeTitleAlignments()
        rsf = ReadSetFilter(0.9)
        self.assertTrue(rsf.accept('title1', titleAlignments))

    def testDuplicateSingleRead(self):
        """
        Testing for acceptance against a read set filter that has already
        seen the exact set should return C{False} if the C{minNew} threshold
        is non-zero.
        """
        rsf = ReadSetFilter(0.9)
        rsf.accept('title1', self.makeTitleAlignments(0))
        self.assertFalse(rsf.accept('title2', self.makeTitleAlignments(0)))

    def testDuplicateSingleReadZeroThreshold(self):
        """
        Testing for acceptance against a read set filter that has already
        seen the exact set should return C{True} if the C{minNew} threshold
        is zero.
        """
        rsf = ReadSetFilter(0.0)
        rsf.accept('title1', self.makeTitleAlignments(0))
        self.assertTrue(rsf.accept('title2', self.makeTitleAlignments(0)))

    def testDifferentSet(self):
        """
        Testing for acceptance against a read set filter that has seen a set
        should return C{True} if the new set is totally different.
        """
        rsf = ReadSetFilter(1.0)
        rsf.accept('title1', self.makeTitleAlignments(0))
        self.assertTrue(rsf.accept('title2', self.makeTitleAlignments(1)))

    def testSufficientlyDifferent(self):
        """
        Testing for acceptance against a read set filter that has seen several
        sets should return C{True} if the new set is sufficiently different.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', self.makeTitleAlignments(0, 1, 2, 3, 4))
        rsf.accept('title2', self.makeTitleAlignments(5, 6, 7, 8, 9))
        self.assertTrue(rsf.accept('title3',
                                   self.makeTitleAlignments(0, 1, 2, 5, 6, 7)))

    def testInsufficientlyDifferent(self):
        """
        Testing for acceptance against a read set filter that has seen several
        sets should return C{False} if the new set is insufficiently different.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', self.makeTitleAlignments(0, 1, 2, 3, 4))
        rsf.accept('title2', self.makeTitleAlignments(5, 6, 7, 8, 9))
        self.assertFalse(rsf.accept('title3',
                                    self.makeTitleAlignments(0, 1, 2, 11)))

    def testThresholdRoundsUp(self):
        """
        Testing for acceptance should round up the needed number of new reads.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', self.makeTitleAlignments(0, 1, 2, 3, 4))
        # If we pass a read set of size three, two of the reads will need to be
        # different.
        self.assertFalse(rsf.accept('title2',
                                    self.makeTitleAlignments(0, 1, 6)))

    def testRepeatTitle(self):
        """
        Testing for acceptance on a title that has been seen before (in an
        accepted read set) must raise C{AssertionError}.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', self.makeTitleAlignments(0, 1, 2, 3, 4))
        self.assertRaises(AssertionError, rsf.accept, 'title1',
                          self.makeTitleAlignments())

    def testInvalidates(self):
        """
        It must be possible to retrieve the list of titles that were
        invalidated by an earlier title's read set.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', self.makeTitleAlignments(0))
        rsf.accept('title2', self.makeTitleAlignments(0))
        rsf.accept('title3', self.makeTitleAlignments(1))
        rsf.accept('title4', self.makeTitleAlignments(0))
        self.assertEqual(['title2', 'title4'], rsf.invalidates('title1'))

    def testInvalidatesEmpty(self):
        """
        The list of titles invalidated by an earlier title that didn't
        invalidate anything must be empty.
        """
        rsf = ReadSetFilter(0.5)
        self.assertEqual([], rsf.invalidates('title1'))


class FakeCursor(object):
    def __init__(self, results):
        self._results = results
        self._index = -1

    def execute(self, p):
        pass

    def fetchone(self):
        self._index += 1
        return self._results[self._index]

    def close(self):
        pass


class FakeDbConnection(object):
    """
    FakeDbConnection and FakeCursor fake results
    for database calls.
    """
    def __init__(self, results):
        self._results = results
        self.open = True

    def cursor(self):
        return FakeCursor(self._results)

    def close(self):
        self.open = False

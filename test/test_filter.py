from unittest import TestCase

from dark.filter import BitScoreFilter, HitInfoFilter, TitleFilter


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


class HitInfoFilterTest(TestCase):
    """
    Tests for the L{dark.filter.HitInfoFilter} class.
    """

    def testNoRestriction(self):
        """
        Testing for acceptance against a hit info filter that has no
        restrictions should return C{True}.
        """
        hif = HitInfoFilter()
        self.assertEqual(True, hif.accept({}))

    def testMinSequenceLen(self):
        """
        Testing for acceptance against a hit info filter with a minSequenceLen
        restriction must work.
        """
        hif = HitInfoFilter(minSequenceLen=10)
        self.assertEqual(True, hif.accept({'length': 12}))
        self.assertEqual(False, hif.accept({'length': 8}))

    def testMaxSequenceLen(self):
        """
        Testing for acceptance against a hit info filter with a maxSequenceLen
        restriction must work.
        """
        hif = HitInfoFilter(maxSequenceLen=100)
        self.assertEqual(True, hif.accept({'length': 90}))
        self.assertEqual(False, hif.accept({'length': 101}))

    def testMinMatchingReads(self):
        """
        Testing for acceptance against a hit info filter with a
        minMatchingReads restriction must work.
        """
        hif = HitInfoFilter(minMatchingReads=5)
        self.assertEqual(True, hif.accept({'readCount': 10}))
        self.assertEqual(False, hif.accept({'readCount': 2}))

    def testMaxMeanEValue(self):
        """
        Testing for acceptance against a hit info filter with a maxMeanEValue
        restriction must work.
        """
        hif = HitInfoFilter(maxMeanEValue=1e-3)
        self.assertEqual(True, hif.accept({'eMean': 1e-4}))
        self.assertEqual(False, hif.accept({'eMean': 1e-2}))

    def testMaxMedianEValue(self):
        """
        Testing for acceptance against a hit info filter with a maxMedianEValue
        restriction must work.
        """
        hif = HitInfoFilter(maxMedianEValue=1e-3)
        self.assertEqual(True, hif.accept({'eMedian': 1e-4}))
        self.assertEqual(False, hif.accept({'eMedian': 1e-2}))

    def testMaxMinEValue(self):
        """
        Testing for acceptance against a hit info filter with a withEBetterThan
        restriction must work.
        """
        hif = HitInfoFilter(withEBetterThan=1e-3)
        self.assertEqual(True, hif.accept({'eMin': 1e-4}))
        self.assertEqual(False, hif.accept({'eMin': 1e-2}))


class BitScoreFilterTest(TestCase):
    """
    Tests for the L{dark.filter.BitScoreFilter} class.
    """

    def testNoRestriction(self):
        """
        Testing for acceptance against a bit score filter that has no
        restrictions should return C{True}.
        """
        bsf = BitScoreFilter()
        self.assertEqual(True, bsf.accept({}))

    def testMinMeanBitScore(self):
        """
        Testing for acceptance against a bit score filter with a
        minMeanBitScore restriction must work.
        """
        bsf = BitScoreFilter(minMeanBitScore=10)
        self.assertEqual(True, bsf.accept({'bitScoreMean': 11}))
        self.assertEqual(False, bsf.accept({'bitScoreMean': 9}))

    def testMinMedianBitScore(self):
        """
        Testing for acceptance against a bit score filter with a
        minMedianBitScore restriction must work.
        """
        bsf = BitScoreFilter(minMedianBitScore=10)
        self.assertEqual(True, bsf.accept({'bitScoreMedian': 11}))
        self.assertEqual(False, bsf.accept({'bitScoreMedian': 9}))

    def testWithBitScoreBetterThan(self):
        """
        Testing for acceptance against a bit score filter with a
        withBitScoreBetterThan restriction must work.
        """
        bsf = BitScoreFilter(withBitScoreBetterThan=10)
        self.assertEqual(True, bsf.accept({'bitScoreMax': 11}))
        self.assertEqual(False, bsf.accept({'bitScoreMax': 9}))

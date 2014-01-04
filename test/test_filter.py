from unittest import TestCase

from dark.filter import HitInfoFilter, TitleFilter


class TitleFilterTest(TestCase):
    """
    Tests for the dark.filter.TitleFilter class.
    """

    def testNoRestriction(self):
        """
        Testing for acceptance against a title filter that has no
        restrictions should return C{True}.
        """
        tf = TitleFilter()
        self.assertEqual(True, tf.accept('hey'))

    def testPositiveRegex(self):
        """
        Testing for acceptance against a title filter with a positive regex
        must work.
        """
        tf = TitleFilter(positiveRegex=r'x+\s')
        self.assertEqual(True, tf.accept('hey xxx you'))
        self.assertEqual(False, tf.accept('hey xxyou'))

    def testNegativeRegex(self):
        """
        Testing for acceptance against a title filter with a negative regex
        must work.
        """
        tf = TitleFilter(negativeRegex=r'x+\s')
        self.assertEqual(False, tf.accept('hey xxx you'))
        self.assertEqual(True, tf.accept('hey xxyou'))

    def testTruncation(self):
        """
        Testing for acceptance against a title filter with title truncation
        in effect must work.
        """
        tf = TitleFilter(truncateAfter=r'virus')
        self.assertEqual(True, tf.accept('polyoma virus one'))
        self.assertEqual(False, tf.accept('polyoma virus two'))

    def testWhitelist(self):
        """
        Testing for acceptance against a title filter with a whitelist
        must work even when a title is ruled out for other violations.
        """
        tf = TitleFilter(whitelist=['always ok'], negativeRegex='ok')
        self.assertEqual(True, tf.accept('always ok'))
        self.assertEqual(False, tf.accept('always ok not'))

    def testBlacklist(self):
        """
        Testing for acceptance against a title filter with a blacklist
        must work.
        """
        tf = TitleFilter(blacklist=['never ok'], positiveRegex='ok')
        self.assertEqual(False, tf.accept('never ok'))

    def testWhitelistTakesPrecedenceOverBlacklist(self):
        """
        Testing for acceptance against a title filter with a whitelist
        and a blacklist that contain the same title must work (the whitelist
        takes precedence).
        """
        tf = TitleFilter(whitelist=['always ok'], blacklist=['always ok'])
        self.assertEqual(True, tf.accept('always ok'))

    def testWhitelistOnly(self):
        """
        Testing for acceptance against a title filter with a whitelist
        and a negative regex that matches everything.
        """
        tf = TitleFilter(whitelist=['always ok'], negativeRegex='.')
        self.assertEqual(True, tf.accept('always ok'))
        self.assertEqual(False, tf.accept('always not ok'))
        self.assertEqual(False, tf.accept('rubbish'))


class HitInfoTest(TestCase):
    """
    Tests for the dark.filter.HitInfoFilter class.
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
        Testing for acceptance against a hit info filter with a maxMinEValue
        restriction must work.
        """
        hif = HitInfoFilter(maxMinEValue=1e-3)
        self.assertEqual(True, hif.accept({'eMin': 1e-4}))
        self.assertEqual(False, hif.accept({'eMin': 1e-2}))

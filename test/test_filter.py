from unittest import TestCase

from dark.filter import HitInfoFilter, ReadSetFilter, TitleFilter


class TitleFilterTest(TestCase):
    """
    Tests for the dark.filter.TitleFilter class.
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
        Testing for acceptance against a hit info filter with a withEBetterThan
        restriction must work.
        """
        hif = HitInfoFilter(withEBetterThan=1e-3)
        self.assertEqual(True, hif.accept({'eMin': 1e-4}))
        self.assertEqual(False, hif.accept({'eMin': 1e-2}))


class ReadSetTest(TestCase):
    """
    Tests for the L{dark.filter.ReadSetFilter} class.
    """

    def testNoRestriction(self):
        """
        Testing for acceptance against a read set filter that has no
        restrictions should return C{True}.
        """
        rsf = ReadSetFilter(0.9)
        self.assertEqual(True, rsf.accept({
            'readNums': set([0])
        }))

    def testDuplicateSingleRead(self):
        """
        Testing for acceptance against a read set filter that has already
        seen the exact set should return C{False} if the C{minNew} threshold
        is non-zero.
        """
        rsf = ReadSetFilter(0.9)
        rsf.accept({
            'readNums': set([0])
        })
        self.assertEqual(False, rsf.accept({
            'readNums': set([0])
        }))

    def testDuplicateSingleReadZeroThreshold(self):
        """
        Testing for acceptance against a read set filter that has already
        seen the exact set should return C{True} if the C{minNew} threshold
        is zero.
        """
        rsf = ReadSetFilter(0.0)
        rsf.accept({
            'readNums': set([0])
        })
        self.assertEqual(True, rsf.accept({
            'readNums': set([0])
        }))

    def testDifferentSet(self):
        """
        Testing for acceptance against a read set filter that has seen a set
        should return C{True} if the new set is totally different.
        """
        rsf = ReadSetFilter(1.0)
        rsf.accept({
            'readNums': set([0])
        })
        self.assertEqual(True, rsf.accept({
            'readNums': set([1])
        }))

    def testSufficientlyDifferent(self):
        """
        Testing for acceptance against a read set filter that has seen several
        sets should return C{True} if the new set is sufficiently different.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept({
            'readNums': set([0, 1, 2, 3, 4])
        })
        rsf.accept({
            'readNums': set([5, 6, 7, 8, 9])
        })
        self.assertEqual(True, rsf.accept({
            'readNums': set([0, 1, 2, 5, 6, 7])
        }))

    def testInsufficientlyDifferent(self):
        """
        Testing for acceptance against a read set filter that has seen several
        sets should return C{False} if the new set is insufficiently different.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept({
            'readNums': set([0, 1, 2, 3, 4])
        })
        rsf.accept({
            'readNums': set([5, 6, 7, 8, 9])
        })
        self.assertEqual(False, rsf.accept({
            'readNums': set([0, 1, 2, 11])
        }))

    def testThresholdRoundsUp(self):
        """
        Testing for acceptance should round up the needed number of new reads.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept({
            'readNums': set([0, 1, 2, 3, 4])
        })
        # If we pass a read set of size three, two of the reads will need to be
        # different.
        self.assertEqual(False, rsf.accept({
            'readNums': set([0, 1, 6])
        }))

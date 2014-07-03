from unittest import TestCase

from dark.filter import (BitScoreFilter, HitInfoFilter, ReadSetFilter,
                         TitleFilter, TaxonomyFilter)
from dark.blast import BlastHits


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
        self.assertEqual(True, rsf.accept('title1', {
            'readNums': set([0])
        }))

    def testDuplicateSingleRead(self):
        """
        Testing for acceptance against a read set filter that has already
        seen the exact set should return C{False} if the C{minNew} threshold
        is non-zero.
        """
        rsf = ReadSetFilter(0.9)
        rsf.accept('title1', {
            'readNums': set([0])
        })
        self.assertEqual(False, rsf.accept('title2', {
            'readNums': set([0])
        }))

    def testDuplicateSingleReadZeroThreshold(self):
        """
        Testing for acceptance against a read set filter that has already
        seen the exact set should return C{True} if the C{minNew} threshold
        is zero.
        """
        rsf = ReadSetFilter(0.0)
        rsf.accept('title1', {
            'readNums': set([0])
        })
        self.assertEqual(True, rsf.accept('title2', {
            'readNums': set([0])
        }))

    def testDifferentSet(self):
        """
        Testing for acceptance against a read set filter that has seen a set
        should return C{True} if the new set is totally different.
        """
        rsf = ReadSetFilter(1.0)
        rsf.accept('title1', {
            'readNums': set([0])
        })
        self.assertEqual(True, rsf.accept('title2', {
            'readNums': set([1])
        }))

    def testSufficientlyDifferent(self):
        """
        Testing for acceptance against a read set filter that has seen several
        sets should return C{True} if the new set is sufficiently different.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', {
            'readNums': set([0, 1, 2, 3, 4])
        })
        rsf.accept('title2', {
            'readNums': set([5, 6, 7, 8, 9])
        })
        self.assertEqual(True, rsf.accept('title3', {
            'readNums': set([0, 1, 2, 5, 6, 7])
        }))

    def testInsufficientlyDifferent(self):
        """
        Testing for acceptance against a read set filter that has seen several
        sets should return C{False} if the new set is insufficiently different.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', {
            'readNums': set([0, 1, 2, 3, 4])
        })
        rsf.accept('title2', {
            'readNums': set([5, 6, 7, 8, 9])
        })
        self.assertEqual(False, rsf.accept('title3', {
            'readNums': set([0, 1, 2, 11])
        }))

    def testThresholdRoundsUp(self):
        """
        Testing for acceptance should round up the needed number of new reads.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', {
            'readNums': set([0, 1, 2, 3, 4])
        })
        # If we pass a read set of size three, two of the reads will need to be
        # different.
        self.assertEqual(False, rsf.accept('title2', {
            'readNums': set([0, 1, 6])
        }))

    def testRepeatTitle(self):
        """
        Testing for acceptance on a title that has been seen before (in an
        accepted read set) must raise C{AssertionError}.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', {
            'readNums': set([0, 1, 2, 3, 4])
        })
        self.assertRaises(AssertionError, rsf.accept, 'title1', {
            'readNums': set([0, 1, 6])
        })

    def testInvalidates(self):
        """
        It must be possible to retrieve the list of titles that were
        invalidated by an earlier title's read set.
        """
        rsf = ReadSetFilter(0.5)
        rsf.accept('title1', {
            'readNums': set([0])
        })
        rsf.accept('title2', {
            'readNums': set([0])
        })
        rsf.accept('title3', {
            'readNums': set([1])
        })
        rsf.accept('title4', {
            'readNums': set([0])
        })
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


class TaxonomyFilterTest(TestCase):
    """
    Tests for L{dark.filter.TaxonomyFilter} class.
    """

    def testNoTaxonomy(self):
        blastHits = BlastHits(None)
        blastHits.addHit('gi|293595919|gb|HM011539.1|', {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage', 'Vira']
        })
        result = blastHits.filterHits(taxonomy=None)
        self.assertEqual(result.titles, {'gi|293595919|gb|HM011539.1|': {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage',
                         'Vira']}})

    def testGivenTaxonomyNotPresent(self):
        blastHits = BlastHits(None)
        blastHits.addHit('gi|293595919|gb|HM011539.1|', {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage',
                         'Vira']
        })
        result = blastHits.filterHits(taxonomy='fiction')
        self.assertEqual(result.titles, {})

    def testGivenTaxonomyPresent(self):
        blastHits = BlastHits(None)
        blastHits.addHit('gi|293595919|gb|HM011539.1|', {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage',
                         'Vira']
        })
        result = blastHits.filterHits(taxonomy='Vira')
        self.assertEqual(result.titles, {'gi|293595919|gb|HM011539.1|': {
                                         'eMin': 0.01,
                                         'taxonomy': ['Merkel cell '
                                                      'polyomavirus',
                                                      'unclassified '
                                                      'Polyomavirus',
                                                      'Polyomavirus',
                                                      'Polyomaviridae',
                                                      'dsDNA viruses, '
                                                      'no RNA stage',
                                                      'Vira']}})

    def test_getTaxIDFromMySql(self):
        gi = 5
        taxonomy = 'Vira'
        db = FakeDbConnection([[15], [2], ['Merkel cell polyomavirus'],
                              [3], ['Polyomavirus'], [2],
                              ['dsDNA viruses'], [1], ['Vira']])
        cursor = db.cursor()
        taxonomyFilter = TaxonomyFilter(gi, taxonomy=taxonomy)
        lineage = taxonomyFilter._getTaxIDFromMySql(cursor, db=db)
        self.assertEqual(['Merkel cell polyomavirus', 'Polyomavirus',
                          'dsDNA viruses', 'Vira'], lineage)

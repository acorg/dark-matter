from unittest import TestCase

from dark.blast import BlastHits


class InterestingRecords(TestCase):
    """
    Tests for the L{blast.BlastHits.filterHits} function.
    """

    def testNoSequences(self):
        """
        If we give C{interesting} no records, we should get nothing back.
        """
        blastHits = BlastHits(None)
        result = blastHits.filterHits()
        self.assertEqual({}, result.titles)

    def testOneSequenceWithNoFiltering(self):
        """
        If we call C{interesting} when we have one title, but provide no
        filtering restrictions we should get the same title back.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'length': 23,
        })
        result = blastHits.filterHits()
        self.assertEqual(
            {
                'a': {
                    'length': 23,
                },
            },
            result.titles)

    def testFilterTitle(self):
        """
        Filtering with a title (i.e., dictionary key) should work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'length': 23,
        })
        blastHits.addHit('b', {
            'length': 15,
        })
        result = blastHits.filterHits(titleRegex='a')
        self.assertEqual(
            {
                'a': {
                    'length': 23,
                },
            },
            result.titles)

    def testFilterTitleCaseInsensitive(self):
        """
        Filtering with a title should work when case doesn't match.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('The ABC Sequence', {
            'length': 23,
        })
        blastHits.addHit('b', {
            'length': 15,
        })
        result = blastHits.filterHits(titleRegex='abc seq')
        self.assertEqual(
            {
                'The ABC Sequence': {
                    'length': 23,
                },
            },
            result.titles)

    def testNegativeFilterTitle(self):
        """
        Filtering negatively with a title (i.e., dictionary key) should work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('The ABC Sequence', {
            'length': 10,
        })
        blastHits.addHit('b', {
            'length': 20,
        })
        result = blastHits.filterHits(negativeTitleRegex='ABC Sequence')
        self.assertEqual(
            {
                'b': {
                    'length': 20,
                }
            },
            result.titles)

    def testNegativeFilterTitleCaseInsensitive(self):
        """
        Filtering negatively with a title (i.e., dictionary key) should
        work even when case doesn't match.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('The ABC Sequence', {
            'length': 10,
        })
        blastHits.addHit('b', {
            'length': 20,
        })
        result = blastHits.filterHits(negativeTitleRegex='abc seq')
        self.assertEqual(
            {
                'b': {
                    'length': 20,
                }
            },
            result.titles)

    def testTruncateTitlesWhenAllMatch(self):
        """
        Truncate titles should give just one result when all sequence titles
        get collapsed.

        Note that 'len' is used in the test because we don't know which of the
        titles will be selected (because iteration order over a dictionary is
        undefined).
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a long virus strain 29744', {
            'length': 23,
        })
        blastHits.addHit('a long virus of unknown aetiology', {
            'length': 33,
        })
        blastHits.addHit('a long virus cryptostrain 9202', {
            'length': 43,
        })
        result = blastHits.filterHits(truncateTitlesAfter='virus')
        self.assertEqual(1, len(result.titles))

    def testTruncateTitlesWhenSomeMatch(self):
        """
        Truncate titles should give the right results when several sequence
        titles are collapsed and some do not.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a long virus strain 29744', {
            'length': 23,
        })
        blastHits.addHit('a long virus of unknown aetiology', {
            'length': 33,
        })
        blastHits.addHit('a long virus cryptostrain 9202', {
            'length': 43,
        })
        blastHits.addHit('something else', {
            'length': 11,
        })
        result = blastHits.filterHits(truncateTitlesAfter='virus')
        self.assertEqual(2, len(result.titles))
        self.assertTrue('something else' in result.titles.keys())

    def testTruncateTitlesWithOneMatch(self):
        """
        Truncate titles shouldn't change anything if only one title matches.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a long virus strain 29744', {
            'length': 23,
        })
        blastHits.addHit('a long virus of unknown aetiology', {
            'length': 33,
        })
        blastHits.addHit('a long virus cryptostrain 9202', {
            'length': 43,
        })
        result = blastHits.filterHits(truncateTitlesAfter='unknown')
        self.assertEqual(
            {
                'a long virus strain 29744': {
                    'length': 23,
                },
                'a long virus of unknown aetiology': {
                    'length': 33,
                },
                'a long virus cryptostrain 9202': {
                    'length': 43,
                }
            },
            result.titles)

    def testTruncateTitlesWithNoMatches(self):
        """
        Truncate titles shouldn't do anything if the titles don't match at all.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a long virus strain 29744', {
            'length': 23,
        })
        blastHits.addHit('a long virus of unknown aetiology', {
            'length': 33,
        })
        blastHits.addHit('a long virus cryptostrain 9202', {
            'length': 43,
        })
        result = blastHits.filterHits(truncateTitlesAfter='xxx')

        self.assertEqual(
            {
                'a long virus strain 29744': {
                    'length': 23,
                },
                'a long virus of unknown aetiology': {
                    'length': 33,
                },
                'a long virus cryptostrain 9202': {
                    'length': 43,
                }
            },
            result.titles)

    def testFilterLengthMinMax(self):
        """
        Filtering with a minimum and maximum sequence length should work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'length': 10
        })
        blastHits.addHit('b', {
            'length': 20
        })
        blastHits.addHit('c', {
            'length': 30
        })
        result = blastHits.filterHits(minSequenceLen=15, maxSequenceLen=25)
        self.assertEqual(
            {
                'b': {
                    'length': 20
                }
            },
            result.titles)

    def testFilterMinMatchingReads(self):
        """
        Filtering with a minimum number of reads should work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'readCount': 10
        })
        blastHits.addHit('b', {
            'readCount': 20
        })
        blastHits.addHit('c', {
            'readCount': 30
        })
        result = blastHits.filterHits(minMatchingReads=15)
        self.assertEqual(
            {
                'b': {
                    'readCount': 20
                },
                'c': {
                    'readCount': 30
                }
            },
            result.titles)

    def testFilterMaxMeanEValue(self):
        """
        Filtering with a maximum mean e-value should work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'eMean': 10
        })
        blastHits.addHit('b', {
            'eMean': 5
        })
        blastHits.addHit('c', {
            'eMean': 30
        })
        result = blastHits.filterHits(maxMeanEValue=6)
        self.assertEqual(
            {
                'b': {
                    'eMean': 5
                }
            },
            result.titles)

    def testFilterMaxMinEValueNoHitsMatch(self):
        """
        Filtering with a maximum minimum e-value should filter out hits
        whose minimal me-values are greater than the passed value.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'eMin': 1e-12
        })
        blastHits.addHit('b', {
            'eMin': 1e-15
        })
        blastHits.addHit('c', {
            'eMin': 1e-5
        })
        result = blastHits.filterHits(withEBetterThan=1e-20)
        self.assertEqual(
            {
            },
            result.titles)

    def testFilterMaxMinEValueAllHitsMatch(self):
        """
        Filtering with a maximum minimum e-value should filter out hits
        whose minimal me-values are greater than the passed value.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'eMin': 1e-12
        })
        blastHits.addHit('b', {
            'eMin': 1e-15
        })
        blastHits.addHit('c', {
            'eMin': 1e-5
        })
        result = blastHits.filterHits(withEBetterThan=1e-3)
        self.assertEqual(
            {
                'a': {
                    'eMin': 1e-12
                },
                'b': {
                    'eMin': 1e-15
                },
                'c': {
                    'eMin': 1e-5
                }
            },
            result.titles)

    def testFilterMaxMinEValueSomeHitsMatch(self):
        """
        Filtering with a maximum minimum e-value should filter out hits
        whose minimal me-values are greater than the passed value.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'eMin': 1e-12
        })
        blastHits.addHit('b', {
            'eMin': 1e-15
        })
        blastHits.addHit('c', {
            'eMin': 1e-5
        })
        result = blastHits.filterHits(withEBetterThan=1e-13)
        self.assertEqual(
            {
                'b': {
                    'eMin': 1e-15
                }
            },
            result.titles)


class TestTitleSorting(TestCase):
    """
    Tests for the L{blast.BlastHits.sortTitles} function.
    """

    def testEmpty(self):
        """
        Sorting when there are no titles must return the empty list.
        """
        blastHits = BlastHits(None)
        result = blastHits.sortTitles('eMean')
        self.assertEqual([], result)

    def testUnknown(self):
        """
        Sorting on an unknown attribute must raise C{ValueError}.
        """
        blastHits = BlastHits(None)
        self.assertRaises(ValueError, blastHits.sortTitles, 'unknown')

    def testEMean(self):
        """
        Sorting on mean e-values must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'eMean': 1e-2,
        })
        blastHits.addHit('ba', {
            'eMean': 1e-4,
        })
        blastHits.addHit('bb', {
            'eMean': 1e-4,
        })
        blastHits.addHit('c', {
            'eMean': 1e-3,
        })
        result = blastHits.sortTitles('eMean')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testEMedian(self):
        """
        Sorting on median e-values must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'eMedian': 1e-2,
        })
        blastHits.addHit('ba', {
            'eMedian': 1e-4,
        })
        blastHits.addHit('bb', {
            'eMedian': 1e-4,
        })
        blastHits.addHit('c', {
            'eMedian': 1e-3,
        })
        result = blastHits.sortTitles('eMedian')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testEMin(self):
        """
        Sorting on minimum e-values must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'eMin': 1e-2,
        })
        blastHits.addHit('ba', {
            'eMin': 1e-4,
        })
        blastHits.addHit('bb', {
            'eMin': 1e-4,
        })
        blastHits.addHit('c', {
            'eMin': 1e-3,
        })
        result = blastHits.sortTitles('eMin')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testBitScoreMean(self):
        """
        Sorting on mean bitscore must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'bitScoreMean': 150,
        })
        blastHits.addHit('ba', {
            'bitScoreMean': 50,
        })
        blastHits.addHit('bb', {
            'bitScoreMean': 50,
        })
        blastHits.addHit('c', {
            'bitScoreMean': 100,
        })
        result = blastHits.sortTitles('bitScoreMean')
        self.assertEqual(['a', 'c', 'ba', 'bb'], result)

    def testBitScoreMax(self):
        """
        Sorting on maximum bitscore must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'bitScoreMax': 150,
        })
        blastHits.addHit('ba', {
            'bitScoreMax': 50,
        })
        blastHits.addHit('bb', {
            'bitScoreMax': 50,
        })
        blastHits.addHit('c', {
            'bitScoreMax': 100,
        })
        result = blastHits.sortTitles('bitScoreMax')
        self.assertEqual(['a', 'c', 'ba', 'bb'], result)

    def testBitScoreMedian(self):
        """
        Sorting on median bitscore must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'bitScoreMedian': 150,
        })
        blastHits.addHit('ba', {
            'bitScoreMedian': 50,
        })
        blastHits.addHit('bb', {
            'bitScoreMedian': 50,
        })
        blastHits.addHit('c', {
            'bitScoreMedian': 100,
        })
        result = blastHits.sortTitles('bitScoreMedian')
        self.assertEqual(['a', 'c', 'ba', 'bb'], result)

    def testReadCount(self):
        """
        Sorting on read count must work, including a secondary sort on title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'readCount': 1,
        })
        blastHits.addHit('ba', {
            'readCount': 5,
        })
        blastHits.addHit('bb', {
            'readCount': 5,
        })
        blastHits.addHit('c', {
            'readCount': 3,
        })
        result = blastHits.sortTitles('readCount')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testLength(self):
        """
        Sorting on sequence length must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'length': 100,
        })
        blastHits.addHit('ba', {
            'length': 500,
        })
        blastHits.addHit('bb', {
            'length': 500,
        })
        blastHits.addHit('c', {
            'length': 300,
        })
        result = blastHits.sortTitles('length')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testTitle(self):
        """
        Sorting on title must work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {})
        blastHits.addHit('ba', {})
        blastHits.addHit('bb', {})
        blastHits.addHit('c', {})
        result = blastHits.sortTitles('title')
        self.assertEqual(['a', 'ba', 'bb', 'c'], result)


class TestTitleSortingOnPlotInfo(TestCase):
    """
    Tests for the L{blast.BlastHits.sortTitlesOnPlotInfo} function.
    """
    def testPlotInfoEmpty(self):
        """
        Sorting when there are no titles must return the empty list.
        """
        blastHits = BlastHits(None)
        result = blastHits.sortTitlesOnPlotInfo('eMean')
        self.assertEqual([], result)

    def testPlotInfoUnknown(self):
        """
        Sorting on an unknown attribute must raise C{ValueError}.
        """
        blastHits = BlastHits(None)
        self.assertRaises(ValueError,
                          blastHits.sortTitlesOnPlotInfo, 'unknown')

    def testPlotInfoEMean(self):
        """
        Sorting on mean e-values must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'plotInfo': {
                'originalEMean': 1e-2,
                }
        })
        blastHits.addHit('ba', {
            'plotInfo': {
                'originalEMean': 1e-4,
                }
        })
        blastHits.addHit('bb', {
            'plotInfo': {
                'originalEMean': 1e-4,
                }
        })
        blastHits.addHit('c', {
            'plotInfo': {
                'originalEMean': 1e-3,
                }
        })
        result = blastHits.sortTitlesOnPlotInfo('eMean')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testPlotInfoEMedian(self):
        """
        Sorting on median e-values must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'plotInfo': {
                'originalEMedian': 1e-2,
                }
        })
        blastHits.addHit('ba', {
            'plotInfo': {
                'originalEMedian': 1e-4,
                }
        })
        blastHits.addHit('bb', {
            'plotInfo': {
                'originalEMedian': 1e-4,
                }
        })
        blastHits.addHit('c', {
            'plotInfo': {
                'originalEMedian': 1e-3,
                }
        })
        result = blastHits.sortTitlesOnPlotInfo('eMedian')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testPlotInfoEMin(self):
        """
        Sorting on minimum e-values must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'plotInfo': {
                'originalEMin': 1e-2,
                }
        })
        blastHits.addHit('ba', {
            'plotInfo': {
                'originalEMin': 1e-4,
                }
        })
        blastHits.addHit('bb', {
            'plotInfo': {
                'originalEMin': 1e-4,
                }
        })
        blastHits.addHit('c', {
            'plotInfo': {
                'originalEMin': 1e-3,
                }
        })
        result = blastHits.sortTitlesOnPlotInfo('eMin')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testPlotInfoBitScoreMax(self):
        """
        Sorting on maximum bitscore must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'plotInfo': {
                'bitScoreMax': 150
                }
        })
        blastHits.addHit('ba', {
            'plotInfo': {
                'bitScoreMax': 50,
                }
        })
        blastHits.addHit('bb', {
            'plotInfo': {
                'bitScoreMax': 50,
                }
        })
        blastHits.addHit('c', {
            'plotInfo': {
                'bitScoreMax': 100,
                }
        })
        result = blastHits.sortTitlesOnPlotInfo('bitScoreMax')
        self.assertEqual(['a', 'c', 'ba', 'bb'], result)

    def testPlotInfoBitScoreMean(self):
        """
        Sorting on mean bitscore must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'plotInfo': {
                'bitScoreMean': 150,
                }
        })
        blastHits.addHit('ba', {
            'plotInfo': {
                'bitScoreMean': 50,
                }
        })
        blastHits.addHit('bb', {
            'plotInfo': {
                'bitScoreMean': 50,
                }
        })
        blastHits.addHit('c', {
            'plotInfo': {
                'bitScoreMean': 100,
                }
        })
        result = blastHits.sortTitlesOnPlotInfo('bitScoreMean')
        self.assertEqual(['a', 'c', 'ba', 'bb'], result)

    def testPlotInfoBitScoreMedian(self):
        """
        Sorting on median bitscore must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'plotInfo': {
                'bitScoreMedian': 150,
                }
        })
        blastHits.addHit('ba', {
            'plotInfo': {
                'bitScoreMedian': 50,
                }
        })
        blastHits.addHit('bb', {
            'plotInfo': {
                'bitScoreMedian': 50,
                }
        })
        blastHits.addHit('c', {
            'plotInfo': {
                'bitScoreMedian': 100,
                }
        })
        result = blastHits.sortTitlesOnPlotInfo('bitScoreMedian')
        self.assertEqual(['a', 'c', 'ba', 'bb'], result)

    def testPlotInfoReadCount(self):
        """
        Sorting on read count must work, including a secondary sort on title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'readCount': 1,
        })
        blastHits.addHit('ba', {
            'readCount': 5,
        })
        blastHits.addHit('bb', {
            'readCount': 5,
        })
        blastHits.addHit('c', {
            'readCount': 3,
        })
        result = blastHits.sortTitlesOnPlotInfo('readCount')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testPlotInfoLength(self):
        """
        Sorting on sequence length must work, including a secondary sort on
        title.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'length': 100,
        })
        blastHits.addHit('ba', {
            'length': 500,
        })
        blastHits.addHit('bb', {
            'length': 500,
        })
        blastHits.addHit('c', {
            'length': 300,
        })
        result = blastHits.sortTitlesOnPlotInfo('length')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testPlotInfoTitle(self):
        """
        Sorting on title must work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {})
        blastHits.addHit('ba', {})
        blastHits.addHit('bb', {})
        blastHits.addHit('c', {})
        result = blastHits.sortTitlesOnPlotInfo('title')
        self.assertEqual(['a', 'ba', 'bb', 'c'], result)

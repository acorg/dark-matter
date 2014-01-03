from unittest import TestCase

from dark.blast import BlastHits


class InterestingRecords(TestCase):
    """
    Tests for the L{blast.BlastHits.interesting} function.
    """

    def testNoSequences(self):
        """
        If we give C{interesting} no records, we should get nothing back.
        """
        blastHits = BlastHits(None)
        result = blastHits.interesting()
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
        result = blastHits.interesting()
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
        result = blastHits.interesting(titleRegex='a')
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
        result = blastHits.interesting(titleRegex='abc seq')
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
        result = blastHits.interesting(negativeTitleRegex='ABC Sequence')
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
        result = blastHits.interesting(negativeTitleRegex='abc seq')
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
        result = blastHits.interesting(truncateTitlesAfter='virus')
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
        result = blastHits.interesting(truncateTitlesAfter='virus')
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
        result = blastHits.interesting(truncateTitlesAfter='unknown')
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
        result = blastHits.interesting(truncateTitlesAfter='xxx')

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
        result = blastHits.interesting(minSequenceLen=15, maxSequenceLen=25)
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
        result = blastHits.interesting(minMatchingReads=15)
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
        result = blastHits.interesting(maxMeanEValue=6)
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
        result = blastHits.interesting(maxMinEValue=1e-20)
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
        result = blastHits.interesting(maxMinEValue=1e-3)
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
        result = blastHits.interesting(maxMinEValue=1e-13)
        self.assertEqual(
            {
                'b': {
                    'eMin': 1e-15
                }
            },
            result.titles)

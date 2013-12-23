from unittest import TestCase

from dark.utils import interestingRecords, filterRecords


class InterestingRecords(TestCase):
    """
    Tests for the utils.interestingRecords function.
    """

    def testNoSequences(self):
        """
        If we give interestingRecords no records, we should get nothing back.
        """
        self.assertEqual({}, interestingRecords({}))

    def testOneSequenceWithNoFiltering(self):
        """
        If we give interestingRecords one record, but no filtering restrictions
        we should get the same record back
        """
        self.assertEqual({'a': 23}, interestingRecords({'a': 23}))

    def testFilterTitle(self):
        """
        Filtering with a title (i.e., dictionary key) should work.
        """
        self.assertEqual(
            {
                'a': 23
            },
            interestingRecords(
                {
                    'a': 23,
                    'b': 15
                },
                titleRegex='a'
            ))

    def testFilterTitleCaseInsensitive(self):
        """
        Filtering with a title (i.e., dictionary key) should work
        even when case doesn't match.
        """
        self.assertEqual(
            {
                'The ABC Sequence': 23
            },
            interestingRecords(
                {
                    'The ABC Sequence': 23,
                    'b': 15
                },
                titleRegex='abc seq'
            ))

    def testNegativeFilterTitle(self):
        """
        Filtering negatively with a title (i.e., dictionary key) should work.
        """
        self.assertEqual(
            {
                'b': {
                    'length': 20,
                }
            },
            interestingRecords(
                {
                    'The ABC Sequence': {
                        'length': 10,
                    },
                    'b': {
                        'length': 20,
                    }
                },
                negativeTitleRegex='ABC Sequence',
                minSequenceLen=5
            ))

    def testNegativeFilterTitleCaseInsensitive(self):
        """
        Filtering negatively with a title (i.e., dictionary key) should
        work even when case doesn't match.
        """
        self.assertEqual(
            {
                'b': {
                    'length': 20,
                }
            },
            interestingRecords(
                {
                    'The ABC Sequence': {
                        'length': 10,
                    },
                    'b': {
                        'length': 20,
                    }
                },
                negativeTitleRegex='abc seq',
                minSequenceLen=5
            ))

    def testFilterLengthMinMax(self):
        """
        Filtering with a minimum and maximum sequence length should work.
        """
        self.assertEqual(
            {
                'b': {
                    'length': 20
                }
            },
            interestingRecords(
                {
                    'a': {
                        'length': 10
                    },
                    'b': {
                        'length': 20
                    },
                    'c': {
                        'length': 30
                    }
                },
                minSequenceLen=15,
                maxSequenceLen=25
            ))

    def testFilterMinMatchingReads(self):
        """
        Filtering with a minimum number of reads should work.
        """
        self.assertEqual(
            {
                'b': {
                    'count': 20
                },
                'c': {
                    'count': 30
                }
            },
            interestingRecords(
                {
                    'a': {
                        'count': 10
                    },
                    'b': {
                        'count': 20
                    },
                    'c': {
                        'count': 30
                    }
                },
                minMatchingReads=15
            ))

    def testFilterMaxMeanEValue(self):
        """
        Filtering with a maximum mean e-value should work.
        """
        self.assertEqual(
            {
                'b': {
                    'eMean': 5
                }
            },
            interestingRecords(
                {
                    'a': {
                        'eMean': 10
                    },
                    'b': {
                        'eMean': 5
                    },
                    'c': {
                        'eMean': 30
                    }
                },
                maxMeanEValue=6
            ))

    def testFilterMaxMinEValueNoHitsMatch(self):
        """
        Filtering with a maximum minimum e-value should filter out hits
        whose minimal me-values are greater than the passed value.
        """
        self.assertEqual(
            {
            },
            interestingRecords(
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
                maxMinEValue=1e-20
            ))

    def testFilterMaxMinEValueAllHitsMatch(self):
        """
        Filtering with a maximum minimum e-value should filter out hits
        whose minimal me-values are greater than the passed value.
        """
        data = {
            'a': {
                'eMin': 1e-12
            },
            'b': {
                'eMin': 1e-15
            },
            'c': {
                'eMin': 1e-5
            }
        }

        self.assertEqual(data, interestingRecords(data, maxMinEValue=1e-3))

    def testFilterMaxMinEValueSomeHitsMatch(self):
        """
        Filtering with a maximum minimum e-value should filter out hits
        whose minimal me-values are greater than the passed value.
        """
        self.assertEqual(
            {
                'b': {
                    'eMin': 1e-15
                }
            },
            interestingRecords(
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
                maxMinEValue=1e-13
            ))


class TestFilterRecords(TestCase):
    """
    Test the filterRecords function.
    """

    def testNoHits(self):
        self.assertEqual({}, filterRecords({}, lambda x: True))

    def testAllTrue(self):
        summary = {
            'a': 234,
            'z': 'hey',
        }
        self.assertEqual(summary, filterRecords(summary, lambda x: True))

    def testAllFalse(self):
        summary = {
            'a': 234,
            'z': 'hey',
        }
        self.assertEqual({}, filterRecords(summary, lambda x: False))

    def testMatchSome1(self):
        summary = {
            'a': {
                'title': 'A sequence',
                'length': 44,
            },
            'b': {
                'title': 'Silly',
                'length': 21
            },
            'c': {
                'title': 'Billy',
                'length': 100
            },
        }
        self.assertEqual(
            {
                'a': {
                    'title': 'A sequence',
                    'length': 44,
                },
            },
            filterRecords(summary, lambda x: x['title'].find('sequence') > -1))

    def testMatchSome2(self):
        summary = {
            'a': {
                'title': 'A sequence',
                'length': 44,
            },
            'b': {
                'title': 'Silly',
                'length': 21
            },
            'c': {
                'title': 'Billy',
                'length': 100
            },
        }
        self.assertEqual(
            {
                'a': {
                    'title': 'A sequence',
                    'length': 44,
                },
                'c': {
                    'title': 'Billy',
                    'length': 100
                },
            },
            filterRecords(summary, lambda x: x['length'] > 30))

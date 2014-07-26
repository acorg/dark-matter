from dark.blast.blast import BlastHits, BlastRecords


class TestFilterHits(object):
    """
    Tests for the L{blast.BlastHits.filterHits} function.
    """

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


class TestTitleSorting(object):
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


class TestTitleSortingOnPlotInfo(object):
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
            'plotInfo': {},
            'readCount': 1,
        })
        blastHits.addHit('ba', {
            'plotInfo': {},
            'readCount': 5,
        })
        blastHits.addHit('bb', {
            'plotInfo': {},
            'readCount': 5,
        })
        blastHits.addHit('c', {
            'plotInfo': {},
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
            'plotInfo': {},
        })
        blastHits.addHit('ba', {
            'length': 500,
            'plotInfo': {},
        })
        blastHits.addHit('bb', {
            'length': 500,
            'plotInfo': {},
        })
        blastHits.addHit('c', {
            'length': 300,
            'plotInfo': {},
        })
        result = blastHits.sortTitlesOnPlotInfo('length')
        self.assertEqual(['ba', 'bb', 'c', 'a'], result)

    def testPlotInfoTitle(self):
        """
        Sorting on title must work.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('a', {
            'plotInfo': {},
        })
        blastHits.addHit('ba', {
            'plotInfo': {},
        })
        blastHits.addHit('bb', {
            'plotInfo': {},
        })
        blastHits.addHit('c', {
            'plotInfo': {},
        })
        result = blastHits.sortTitlesOnPlotInfo('title')
        self.assertEqual(['a', 'ba', 'bb', 'c'], result)


class TestComputePlotInfo(object):
    """
    Tests for the L{blast.BlastHits.computePlotInfo} function.
    """

    def testSortAfterComputePlotInfo(self):
        """
        If there is no plot info available for some titles after calling
        computePlotInfo, we should still be able to call the various
        plotInfo-based sort routines.
        """
        params = {
            'application': 'BLASTN',
        }
        record1 = {
            'query': 'a',
            'alignments': [
                {
                    'length': 37108,
                    'hsps': [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 1.0,
                            'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    'title': 'title-a',
                }
            ]
        }

        record2 = {
            'query': 'b',
            'alignments': [
                {
                    'length': 37108,
                    'hsps': [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 3.0,
                            'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    'title': 'title-b',
                }
            ]
        }
        mockOpener = mockOpen(read_data=dumps(params) + '\n' +
                              dumps(record1) + '\n' +
                              dumps(record2) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('a.json')
            blastHits = BlastHits(blastRecords)
            blastHits.addHit('title-a', {
                'length': 100,
                'readNums': [0],
            })
            blastHits.addHit('title-b', {
                'length': 100,
                'readNums': [1],
            })
            # Fake out the reading of the FASTA file.
            blastHits.fasta = ['a' * 68, 'a' * 68]
            blastHits.computePlotInfo()
            # And we must be able to sort without error.
            self.assertEqual(['title-a', 'title-b'],
                             blastHits.sortTitlesOnPlotInfo('eMin'))


class TestTaxonomyFiltering:
    def testAllTaxonomy(self):
        blastHits = BlastHits(None)
        blastHits.addHit('gi|293595919|gb|HM011539.1|', {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage', 'Vira']
        })
        result = blastHits.filterHits(taxonomy='all')
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

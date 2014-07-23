from dark.blast.blast import BlastHits, BlastRecords


class TestBlastRecords(object):
    """
    Test the reading of BLAST records.
    """

    def testECutoffRemovesEntireAlignment(self):
        """
        If the L{BlastRecords} records function is supposed to filter on
        eCutoff and the cut-off value results in an alignment with no HSPs,
        then the alignment must be removed entirely.
        """
        params = {
            'application': 'BLASTN',
        }
        record1 = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25854e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 182.092,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 2220,
                    "hsps": [
                        {
                            "sbjct_end": 1841,
                            "expect": 1.25854e-43,
                            "sbjct": "AATCCAGGGAATCTAATAAAATAATCAA",
                            "sbjct_start": 1942,
                            "query": "AATCCAGGGAATCTTAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 180,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(params) + '\n' +
                              dumps(record1) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            records = list(BlastRecords('file.json').records(
                eCutoff=1e-15))
            self.assertEqual(1, len(records))
            self.assertEqual(1, len(records[0].alignments))
            self.assertEqual('Merkel2', records[0].alignments[0].title)

    def testECutoffRemovesHsps(self):
        """
        If the L{BlastRecords} records function is supposed to filter on
        eCutoff and the cut-off value results in some HSPs being invalid,
        then those HSPs must be removed entirely.
        """
        params = {
            'application': 'BLASTN',
        }
        record1 = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25854e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 182.092,
                            "query_start": 362
                        },
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25e-20,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 182.092,
                            "query_start": 0
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 2220,
                    "hsps": [
                        {
                            "sbjct_end": 1841,
                            "expect": 1.25e-43,
                            "sbjct": "AATCCAGGGAATCTAATAAAATAATCAA",
                            "sbjct_start": 1942,
                            "query": "AATCCAGGGAATCTTAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 180,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(params) + '\n' +
                              dumps(record1) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            records = list(BlastRecords('file.json').records(
                eCutoff=1e-15))
            self.assertEqual(1, len(records))
            self.assertEqual(2, len(records[0].alignments))
            # There should only be one HSP left in the alignments for the
            # first read, and it should have the right expect value.
            self.assertEqual(1, len(records[0].alignments[0].hsps))
            self.assertEqual(1.25e-20, records[0].alignments[0].hsps[0].expect)
            # The second alignment should also be present.
            self.assertEqual(1, len(records[0].alignments[1].hsps))
            self.assertEqual(1.25e-43, records[0].alignments[1].hsps[0].expect)

    def testXMLInput(self):
        """
        Test conversion of a chunk of BLAST XML.  This is completely
        minimal, it could assert dozens of other things. But it's not
        really needed.
        """
        record = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" \
"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.2.27+</BlastOutput_version>
  <BlastOutput_reference>[Snip]</BlastOutput_reference>
  <BlastOutput_db>virus-nt-20130719</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>ICUR3MX01AGWKS</BlastOutput_query-def>
  <BlastOutput_query-len>63</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>2</Parameters_sc-match>
      <Parameters_sc-mismatch>-3</Parameters_sc-mismatch>
      <Parameters_gap-open>5</Parameters_gap-open>
      <Parameters_gap-extend>2</Parameters_gap-extend>
      <Parameters_filter>L;m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>ICUR3MX01AGWKS</Iteration_query-def>
      <Iteration_query-len>63</Iteration_query-len>
      <Iteration_hits></Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1331289</Statistics_db-num>
          <Statistics_db-len>1634236772</Statistics_db-len>
          <Statistics_hsp-len>28</Statistics_hsp-len>
          <Statistics_eff-space>55893623800</Statistics_eff-space>
          <Statistics_kappa>0.41</Statistics_kappa>
          <Statistics_lambda>0.625</Statistics_lambda>
          <Statistics_entropy>0.78</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
      <Iteration_message>No hits found</Iteration_message>
    </Iteration>
    <Iteration>
      <Iteration_iter-num>2</Iteration_iter-num>
      <Iteration_query-ID>Query_2</Iteration_query-ID>
      <Iteration_query-def>ICUR3MX01C58U5</Iteration_query-def>
      <Iteration_query-len>21</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|157694982|gb|EF990688.1|</Hit_id>
          <Hit_def>Rotavirus A strain RVA/Pig-tc/VEN/A131/1988</Hit_def>
          <Hit_accession>EF990688</Hit_accession>
          <Hit_len>1059</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>30.1402</Hsp_bit-score>
              <Hsp_score>32</Hsp_score>
              <Hsp_evalue>4.0824</Hsp_evalue>
              <Hsp_query-from>6</Hsp_query-from>
              <Hsp_query-to>21</Hsp_query-to>
              <Hsp_hit-from>738</Hsp_hit-from>
              <Hsp_hit-to>753</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>16</Hsp_identity>
              <Hsp_positive>16</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>16</Hsp_align-len>
              <Hsp_qseq>ATTCATCAGTAGCAAT</Hsp_qseq>
              <Hsp_hseq>ATTCATCAGTAGCAAT</Hsp_hseq>
              <Hsp_midline>||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gi|157694972|gb|EF990692.1|</Hit_id>
          <Hit_def>Rotavirus A strain RVA/Pig-tc/VEN/A411/1989</Hit_def>
          <Hit_accession>EF990692</Hit_accession>
          <Hit_len>1059</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>30.1402</Hsp_bit-score>
              <Hsp_score>32</Hsp_score>
              <Hsp_evalue>4.0824</Hsp_evalue>
              <Hsp_query-from>6</Hsp_query-from>
              <Hsp_query-to>21</Hsp_query-to>
              <Hsp_hit-from>738</Hsp_hit-from>
              <Hsp_hit-to>753</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>16</Hsp_identity>
              <Hsp_positive>16</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>16</Hsp_align-len>
              <Hsp_qseq>ATTCATCAGTAGCAAT</Hsp_qseq>
              <Hsp_hseq>ATTCATCAGTAGCAAT</Hsp_hseq>
              <Hsp_midline>||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1331289</Statistics_db-num>
          <Statistics_db-len>1634236772</Statistics_db-len>
          <Statistics_hsp-len>18</Statistics_hsp-len>
          <Statistics_eff-space>4830820710</Statistics_eff-space>
          <Statistics_kappa>0.41</Statistics_kappa>
          <Statistics_lambda>0.625</Statistics_lambda>
          <Statistics_entropy>0.78</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
        """
        mockOpener = mockOpen(read_data=record)
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('file.xml', None, None)
            record1, record2 = list(blastRecords.records())
            self.assertEqual(0, len(record1.alignments))
            self.assertEqual(2, len(record2.alignments))


class TestFilterHits(object):
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

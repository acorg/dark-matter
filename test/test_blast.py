from json import dumps
from mock import patch
from unittest import TestCase

from dark.blast import BlastHits, BlastRecords
from mocking import mockOpen


class TestBlastRecords(TestCase):
    """
    Test the reading of BLAST records.

    TODO: This class is highly incomplete.
    """

    def testEmptyJSONInput(self):
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('file.json', None, None)
            self.assertEqual([], list(blastRecords.records()))

    def testNonJSONInput(self):
        """
        When given a file whose contents are not JSON, attempting to
        read the BLAST records from it must raise a C{ValueError}.
        """
        mockOpener = mockOpen(read_data='not JSON')
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('file.json', None, None)
            self.assertRaises(ValueError, list, blastRecords.records())

    def testEmptyXMLInput(self):
        """
        When an XML input file is empty, BioPython's C{NCBIXML.parse}
        raises a C{ValueError}.
        """
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('file.xml', None, None)
            generator = blastRecords.records()
            self.assertRaises(ValueError, list, generator)

    def testJSONParams(self):
        """
        BLAST parameters must be found in an input JSON file and stored
        into the C{BlastRecords} instance.
        """
        params = {
            'application': 'BLASTN',
            'blast_cutoff': [
                None,
                None
            ],
            'database': 'manx-shearwater',
            'database_length': 17465129,
            'database_letters': None,
            'database_name': [],
            'database_sequences': 70016,
            'date': '',
            'dropoff_1st_pass': [
                None,
                None
            ],
            'effective_database_length': None,
            'effective_hsp_length': 22,
            'effective_query_length': None,
            'effective_search_space': 382194648.0,
            'effective_search_space_used': None,
            'frameshift': [
                None,
                None
            ],
            'gap_penalties': [
                5,
                2
            ],
            'gap_trigger': [
                None,
                None
            ],
            'gap_x_dropoff': [
                None,
                None
            ],
            'gap_x_dropoff_final': [
                None,
                None
            ],
            'gapped': 0,
            'hsps_gapped': None,
            'hsps_no_gap': None,
            'hsps_prelim_gapped': None,
            'hsps_prelim_gapped_attemped': None,
            'ka_params': [
                0.625,
                0.41,
                0.78
            ],
            'ka_params_gap': [
                None,
                None,
                None
            ],
            'matrix': '',
            'num_good_extends': None,
            'num_hits': None,
            'num_letters_in_database': 17465129,
            'num_seqs_better_e': None,
            'num_sequences': None,
            'num_sequences_in_database': 70016,
            'posted_date': [],
            'query': 'GZG3DGY01ASHXW',
            'query_id': 'Query_1',
            'query_length': 46,
            'query_letters': 46,
            'reference': 'Stephen F. Altschul, Thomas L. Madden, ...',
            'sc_match': 2,
            'sc_mismatch': -3,
            'threshold': None,
            'version': '2.2.28+',
            'window_size': None
        }
        mockOpener = mockOpen(read_data=dumps(params) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('file.json', None, None)
            records = blastRecords.records()
            # There are no records, so calling records.next will raise.
            self.assertRaises(StopIteration, records.next)
            self.assertEqual(params, blastRecords.blastParams)

    def testOneJSONInput(self):
        """
        If a JSON file contains a parameters section and one record, it must be
        read correctly.
        """
        params = {
            'application': 'BLASTN',
        }

        record = {
            'query': 'ICUR3MX01C6VDK',
            'alignments': [
                {
                    'length': 37108,
                    'hsps': [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 3.29804,
                            'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    'title': 'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                }
            ]
        }
        mockOpener = mockOpen(
            read_data=dumps(params) + '\n' + dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('file.json', None, None)
            self.assertEqual(1, len(list(blastRecords.records())))

    def testLimitZero(self):
        """
        If the L{BlastRecords} is limited to reading zero records that limit
        must be respected.
        """
        params = {
            'application': 'BLASTN',
        }

        record = {
            'query': 'a',
            'alignments': [
                {
                    'length': 37108,
                    'hsps': [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 3.29804,
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
        mockOpener = mockOpen(
            read_data=dumps(params) + '\n' + dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            records = list(BlastRecords('file.json', limit=0).records())
            self.assertEqual(0, len(records))

    def testLimitOne(self):
        """
        If the L{BlastRecords} is limited to reading a certain number of
        records that limit must be respected.
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
                            'expect': 3.29804,
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
                            'expect': 3.29804,
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
            records = list(BlastRecords('file.json', limit=1).records())
            self.assertEqual(1, len(records))
            self.assertEqual('title-a', records[0].descriptions[0].title)

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


class TestFilterHits(TestCase):
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

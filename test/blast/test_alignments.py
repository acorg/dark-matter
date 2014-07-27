# TODO: Add tests based on taxonomy, once we know how to mock mysql.

from copy import deepcopy
import bz2
from json import dumps
from unittest import TestCase
from ..mocking import mockOpen
from mock import patch

from dark.reads import Read, Reads
from dark.alignments import Alignment, ReadAlignments, bestAlignment
from dark.blast.hsp import EValueHSP
from dark.blast.alignments import (
    BlastReadsAlignments, numericallySortFilenames)


# Sample BLAST parameters.
PARAMS = {
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

RECORD0 = {
    'query': 'id0',
    'alignments': [
        {
            'length': 37000,
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
        },
        {
            'length': 38000,
            'hsps': [
                {
                    'bits': 20,
                    'sbjct_end': 12400,
                    'expect': 3.29804,
                    'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                    'sbjct_start': 12362,
                    'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                    'frame': [1, 1],
                    'query_end': 68,
                    'query_start': 28
                }
            ],
            'title': 'gi|887699|gb|DQ37780 Squirrelpox virus 55',
        }
    ]
}

RECORD1 = {
    'query': 'id1',
    'alignments': [
        {
            'length': 35000,
            'hsps': [
                {
                    'bits': 20,
                    'sbjct_end': 11400,
                    'expect': 3.29804,
                    'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                    'sbjct_start': 11362,
                    'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                    'frame': [1, 1],
                    'query_end': 68,
                    'query_start': 28
                }
            ],
            'title': 'gi|887699|gb|DQ37780 Monkeypox virus 456',
        },
        {
            'length': 35000,
            'hsps': [
                {
                    'bits': 20,
                    'sbjct_end': 10400,
                    'expect': 3.29804,
                    'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                    'sbjct_start': 10362,
                    'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                    'frame': [1, 1],
                    'query_end': 68,
                    'query_start': 28
                }
            ],
            'title': 'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',
        }
    ]
}

RECORD2 = {
    'query': 'id2',
    'alignments': [
        {
            'length': 30000,
            'hsps': [
                {
                    'bits': 20,
                    'sbjct_end': 1400,
                    'expect': 3.29804,
                    'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                    'sbjct_start': 1362,
                    'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                    'frame': [1, 1],
                    'query_end': 68,
                    'query_start': 28
                }
            ],
            'title': 'gi|887699|gb|DQ37780 Cowpox virus 15',
        }
    ]
}


class BZ2(object):
    """
    A BZ2File mock.
    """
    def __init__(self, data):
        self._data = data
        self._index = 0

    def close(self):
        pass

    def readline(self):
        self._index += 1
        return self._data[self._index - 1]

    def __iter__(self):
        index = self._index
        self._index = len(self._data)
        return iter(self._data[index:])


class TestNumericallySortFilenames(TestCase):
    """
    Test the numericallySortFilenames function.
    """

    def testNoNames(self):
        """
        An empty list must be returned when an empty list is given.
        """
        self.assertEqual([], numericallySortFilenames([]))

    def testOneNonNumericName(self):
        """
        A list with a single non-numeric name should result in that same
        name being returned.
        """
        self.assertEqual(['hey'], numericallySortFilenames(['hey']))

    def testOneNumericName(self):
        """
        A list with a single numeric name should result in that same
        name being returned.
        """
        self.assertEqual(['3.json'], numericallySortFilenames(['3.json']))

    def testSeveralNames(self):
        """
        A list with several numeric names should result in a correctly
        sorted list of names being returned.
        """
        self.assertEqual(
            ['1.json', '2.json', '3.json'],
            numericallySortFilenames(['3.json', '1.json', '2.json']))

    def testSeveralNamesWithUnequalPrefixLengths(self):
        """
        A list with several numeric names whose numeric prefixes differ
        in length should result in a correctly sorted list of names being
        returned.
        """
        self.assertEqual(
            ['2.json', '3.json', '21.json', '35.json', '250.json'],
            numericallySortFilenames(
                ['3.json', '21.json', '35.json', '250.json', '2.json']))


class TestBlastReadsAlignments(TestCase):
    """
    Test the BlastReadsAlignments class.
    """

    def testEmptyJSONInput(self):
        """
        When a JSON input file is empty, a C{ValueError} must be raised
        on trying to read it.
        """
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            error = "JSON file 'file.json' was empty."
            self.assertRaisesRegexp(
                ValueError, error, BlastReadsAlignments, reads, 'file.json')

    def testNonJSONInput(self):
        """
        When given a file whose contents are not JSON, attempting to
        read the BLAST hits from it must raise a C{ValueError}.
        """
        mockOpener = mockOpen(read_data='not JSON\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            error = ("Could not convert first line of 'file.json' to JSON "
                     "\(No JSON object could be decoded\). "
                     "Line is 'not JSON'.")
            self.assertRaisesRegexp(
                ValueError, error, BlastReadsAlignments, reads, 'file.json')

    def testJSONParams(self):
        """
        BLAST parameters must be extracted from the input JSON file and stored
        correctly.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual(PARAMS, readsAlignments.params)

    def testJSONParamsButNoHits(self):
        """
        BLAST parameters are present but there are no hits, the __iter__
        method of a L{BlastReadsAlignments} instance must not yield anything.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual([], list(readsAlignments))

    def testNotEnoughReads(self):
        """
        If a JSON file contains a parameters section and one hit, but there
        is no read to go with the hit, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(
            read_data=dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            error = ("Read generator failed to yield read number 1 during "
                     "parsing of BLAST file 'file\.json'\.")
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertRaisesRegexp(ValueError, error, list, readsAlignments)

    def testTooManyReads(self):
        """
        If a JSON file contains a parameters section and one hit, but there
        is more than one read, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(
            read_data=dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'G' * 70))
            error = ("Reads iterator contained more reads than the number of "
                     "BLAST records found \(1\)\. First unknown read id is "
                     "'id1'\.")
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertRaisesRegexp(ValueError, error, list, readsAlignments)

    def testIncorrectReadId(self):
        """
        If the query id of a hit does not match the id of the corresponding
        input read, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(
            read_data=dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('not id0', 'A' * 70))
            error = ("The reads you have provided do not match the BLAST "
                     "output: BLAST record query id \(id0\) does "
                     "not match the id of the supposedly corresponding read "
                     "\(not id0\)\.")
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertRaisesRegexp(ValueError, error, list, readsAlignments)

    def testOneCompressedJSONInput(self):
        """
        If a compressed (bz2) JSON file contains a parameters section and one
        record, it must be read correctly.
        """
        result = BZ2([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.return_value = result
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json.bz2')
            self.assertEqual(1, len(list(readsAlignments)))

    def testTwoCompressedJSONInputs(self):
        """
        If two compressed (bz2) JSON files are passed to
        L{BlastReadsAlignments} each with a parameters section and one
        record, both records must be read correctly and the result should
        have 2 records.
        """

        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])
                else:
                    return BZ2([dumps(PARAMS) + '\n', dumps(RECORD1) + '\n'])

        sideEffect = SideEffect()
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, ['file1.json.bz2', 'file2.json.bz2'])
            result = list(readsAlignments)
            self.assertEqual(2, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)

    def testThreeCompressedJSONInputs(self):
        """
        If three compressed (bz2) JSON files are passed to
        L{BlastReadsAlignments} with names that have a numeric prefix and
        each with a parameters section and one record, all records must be
        read correctly and the result should have 3 records in the correct
        order.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename):
                if self.count == 0:
                    self.test.assertEqual('1.json.bz2', filename)
                    self.count += 1
                    return BZ2([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])
                elif self.count == 1:
                    self.test.assertEqual('2.json.bz2', filename)
                    self.count += 1
                    return BZ2([dumps(PARAMS) + '\n', dumps(RECORD1) + '\n'])
                else:
                    self.test.assertEqual('3.json.bz2', filename)
                    return BZ2([dumps(PARAMS) + '\n', dumps(RECORD2) + '\n'])

        sideEffect = SideEffect(self)
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))

            # Note the files are given out of order. Their names will be
            # sorted before they are opened. The sorting of the names is
            # verified in the SideEffect class, above.
            readsAlignments = BlastReadsAlignments(
                reads, ['3.json.bz2', '1.json.bz2', '2.json.bz2'])
            result = list(readsAlignments)
            self.assertEqual(3, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)
            self.assertEqual('id2', result[2].read.id)

    def testIncompatibleParameters(self):
        """
        If two compressed (bz2) JSON files with incompatible parameters
        are given to L{BlastReadsAlignments}, a C{ValueError} must be
        raised when the files are read.
        """

        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])
                else:
                    params = deepcopy(PARAMS)
                    params['application'] = 'Skype'
                    return BZ2([dumps(params) + '\n', dumps(RECORD1) + '\n'])

        sideEffect = SideEffect()
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            error = ("Incompatible BLAST parameters found. The parameters in "
                     "file2.json.bz2 differ from those originally found in "
                     "file1.json.bz2. Summary of differences:\n\tParam "
                     "u'application' initial value u'BLASTN' differs from "
                     "later value u'Skype'")
            readsAlignments = BlastReadsAlignments(
                reads, ['file1.json.bz2', 'file2.json.bz2'])
            self.assertRaisesRegexp(ValueError, error, list, readsAlignments)


class TestBlastReadsAlignmentsFiltering(TestCase):
    """
    Test the BlastReadsAlignments class filter function.
    """

    def testNoResultNoFilteringArgs(self):
        """
        If the L{BlastReadsAlignments} filter function is called with no
        arguments, and there are no hits, it should produce a generator
        that yields no result.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter())
            self.assertEqual(0, len(result))

    def testOneHitNoFilteringArgs(self):
        """
        If the L{BlastReadsAlignments} filter function is called with no
        arguments, and there is one hit, it should produce a generator that
        yields that hit.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter())
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)

    def testLimitZero(self):
        """
        If L{BlastReadsAlignments} is limited to zero result, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(limit=0))
            self.assertEqual(0, len(result))

    def testLimitOne(self):
        """
        If L{BlastReadsAlignments} is limited to one hit, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n' +
                              dumps(RECORD1) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(limit=1))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)

    def testOneAlignmentPerRead(self):
        """
        If L{BlastReadsAlignments} is asked to deliver only the best alignment
        for each read, that must be respected.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25854e-43,
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
                    "title":"Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(oneAlignmentPerRead=True))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual('Merkel1', result[0].alignments[0].hitTitle)

    def testScoreCutoffRemovesEntireAlignment_Bits(self):
        """
        If the L{BlastReadsAlignments} filter function is supposed to filter on
        a scoreCutoff (bit score) and the cut-off value results in an
        alignment with no HSPs, then the alignment must be removed entirely.
        """
        record = {
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
                            "bits": 150,
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

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(scoreCutoff=160))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual('Merkel2', result[0].alignments[0].hitTitle)

    def testScoreCutoffRemovesEntireAlignment_EValue(self):
        """
        If the L{BlastReadsAlignments} filter function is supposed to filter on
        a scoreCutoff (bit score) and the cut-off value results in an
        alignment with no HSPs, then the alignment must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 150,
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
                            "expect": 1.25e-30,
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

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreType='e values')
            result = list(readsAlignments.filter(scoreCutoff=1e-20))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual('Merkel2', result[0].alignments[0].hitTitle)

    def testScoreCutoffRemovesHsps_Bits(self):
        """
        If the L{BlastRecords} records function is supposed to filter on
        scoreCutoff (bit score) and the cut-off value results in some HSPs
        being invalid, then those HSPs must be removed entirely.
        """
        record = {
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
                            "bits": 150,
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
                            "bits": 170,
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

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(scoreCutoff=160))

            # There should only be one HSP left in the alignments for the
            # first read, and it should have the right score (bit score).
            self.assertEqual(1, len(result[0].alignments[0].hsps))
            self.assertEqual(170, result[0].alignments[0].hsps[0].score)

            # The second alignment should also be present.
            self.assertEqual(1, len(result[0].alignments[1].hsps))
            self.assertEqual(180, result[0].alignments[1].hsps[0].score)

    def testScoreCutoffRemovesHsps_EValue(self):
        """
        If the L{BlastRecords} records function is supposed to filter on
        scoreCutoff (bit score) and the cut-off value results in some HSPs
        being invalid, then those HSPs must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 150,
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
                            "bits": 170,
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
                            "expect": 1.25e-30,
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

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreType='e values')
            result = list(readsAlignments.filter(scoreCutoff=1e-15))

            # There should only be one HSP left in the alignments for the
            # first read, and it should have the right score (e-value).
            self.assertEqual(1, len(result[0].alignments[0].hsps))
            self.assertEqual(1.25e-20, result[0].alignments[0].hsps[0].score)

            # The second alignment should also be present.
            self.assertEqual(1, len(result[0].alignments[1].hsps))
            self.assertEqual(1.25e-30, result[0].alignments[1].hsps[0].score)

    def testTitleByRegexCaseInvariant(self):
        """
        Filtering with a title regex must work independent of case.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(titleRegex='sqUIRRel'))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0].alignments[0].hitTitle)

    def testTitleByRegexAllAlignments(self):
        """
        Filtering with a title regex must work in the case that all alignments
        for a hit match the regex.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(titleRegex='squirrel'))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0].alignments[0].hitTitle)

    def testTitleByRegexOneAlignments(self):
        """
        Filtering with a title regex must work in the case that only some
        alignments for a hit match the regex.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(titleRegex='Mummy'))
            self.assertEqual(1, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',
                             result[0].alignments[0].hitTitle)

    def testTitleByNegativeRegexAllAlignments(self):
        """
        Filtering with a negative title regex must work in the case that
        all alignments for a hit are ruled out (in which case the whole hit,
        for read id0 hitting Squirrelpox in this case, must be skipped).
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(
                negativeTitleRegex='squirrel'))
            self.assertEqual(2, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Monkeypox virus 456',
                             result[0].alignments[0].hitTitle)
            self.assertEqual('id2', result[1].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Cowpox virus 15',
                             result[1].alignments[0].hitTitle)

    def testTitleByNegativeRegexOneAlignment(self):
        """
        Filtering with a negative title regex must work in the case that only
        some alignments for a hit are ruled out (in which case only those
        alignments must be removed but the hit is still valid).
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(negativeTitleRegex='Mummy'))
            self.assertEqual(3, len(result))
            self.assertEqual('id1', result[1].read.id)
            self.assertEqual(1, len(result[1].alignments))
            self.assertEqual('gi|887699|gb|DQ37780 Monkeypox virus 456',
                             result[1].alignments[0].hitTitle)

    def testTitleByNegativeRegexOneAlignment(self):
        """
        Filtering with a negative title regex must work in the case that only
        some alignments for a hit are ruled out (in which case only those
        alignments must be removed but the hit is still valid).
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(negativeTitleRegex='Mummy'))
            self.assertEqual(3, len(result))
            self.assertEqual('id1', result[1].read.id)
            self.assertEqual(1, len(result[1].alignments))
            self.assertEqual('gi|887699|gb|DQ37780 Monkeypox virus 456',
                             result[1].alignments[0].hitTitle)

    def testTitleByNegativeRegexMatchesAll(self):
        """
        Filtering with a negative title regex that matches all alignments
        must remove everything and return an empty result.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(negativeTitleRegex='pox'))
            self.assertEqual(0, len(result))

    def testTitleByNegativeRegexMatchingAllWithWhitelist(self):
        """
        Filtering with a negative title regex that matches all alignments
        must remove everything and result in no hits, except for any
        whitelisted titles.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            title = 'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99'
            result = list(readsAlignments.filter(negativeTitleRegex='pox',
                                                 whitelist=[title]))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual(title, result[0].alignments[0].hitTitle)

    def testTitleByRegexMatchingAllWithBlacklist(self):
        """
        Filtering with a title regex that matches all alignments
        must keep everything, except for any blacklisted titles.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            blacklist = ['gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                         'gi|887699|gb|DQ37780 Squirrelpox virus 55']
            result = list(readsAlignments.filter(titleRegex='pox',
                                                 blacklist=blacklist))
            self.assertEqual(2, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual('id2', result[1].read.id)

    def testTitleTruncation(self):
        """
        When truncating titles, if a set of matched sequences has titles that
        are identical up to the truncation word, only the first found is
        returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(truncateTitlesAfter='virus'))
            self.assertEqual(3, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0].alignments))
            # The Squirrelpox virus 55 hit in RECORD0 is not returned.
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0].alignments[0].hitTitle)

    def testMinTitleSequenceLength(self):
        """
        It must be possible to filter alignments based on minimum hit sequence
        length.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(minSequenceLen=37500))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 55',
                             result[0].alignments[0].hitTitle)

    def testMinTitleSequenceLengthNoHits(self):
        """
        It must be possible to filter alignments based on minimum hit sequence
        length and if nothing sufficiently long matches, an empty list of
        alignments must be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(minSequenceLen=1000000))
            self.assertEqual(0, len(result))

    def testMaxTitleSequenceLength(self):
        """
        It must be possible to filter alignments based on maximum hit sequence
        length.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(maxSequenceLen=31000))
            self.assertEqual(1, len(result))
            self.assertEqual('id2', result[0].read.id)
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual('gi|887699|gb|DQ37780 Cowpox virus 15',
                             result[0].alignments[0].hitTitle)

    def testMaxTitleSequenceLengthNoHits(self):
        """
        It must be possible to filter alignments based on maximum hit sequence
        length and if no sufficiently short sequences match, an empty
        list of alignments must be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(maxSequenceLen=10000))
            self.assertEqual(0, len(result))

    def testMinAndMaxTitleSequenceLength(self):
        """
        It must be possible to filter alignments simultaneously on minimum and
        maximum hit sequence length.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(minSequenceLen=37000,
                                                 maxSequenceLen=38000))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(2, len(result[0].alignments))
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0].alignments[0].hitTitle)
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 55',
                             result[0].alignments[1].hitTitle)

    def testMinStart(self):
        """
        It must be possible to filter alignments based on minimum offset in
        the hit sequence.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(minStart=15300))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0].alignments[0].hitTitle)

    def testMinStartNoHits(self):
        """
        It must be possible to filter alignments based on minimum offset in
        the hit sequence, and if no hsps match then an empty result set
        must be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(minStart=100000))
            self.assertEqual(0, len(result))

    def testMaxStop(self):
        """
        It must be possible to filter alignments based on maximum offset in
        the hit sequence.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(maxStop=1500))
            self.assertEqual(1, len(result))
            self.assertEqual('id2', result[0].read.id)
            self.assertEqual(1, len(result[0].alignments))
            self.assertEqual('gi|887699|gb|DQ37780 Cowpox virus 15',
                             result[0].alignments[0].hitTitle)

    def testMaxStopNoHits(self):
        """
        It must be possible to filter alignments based on maximum offset in
        the hit sequence, and if no hsps match then an empty result set must
        be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(maxStop=100))
            self.assertEqual(0, len(result))

    def testMinStartAndMaxstop(self):
        """
        It must be possible to filter alignments based simultaneously on
        mininum and maximum offset in the hit sequence.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(minStart=9000, maxStop=12000))
            self.assertEqual(1, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual(2, len(result[0].alignments))


class TestBestAlignment(TestCase):
    """
    Test the bestAlignment function using EValueHSPs (i.e., where lower scores
    are better).
    """

    def testOneAlignment(self):
        """
        When one alignment is present that alignment must be returned by
        bestAlignment.
        """
        alignment = Alignment(44, 'Seq 1')
        alignment.addHsp(EValueHSP(score=10))
        alignment.addHsp(EValueHSP(score=9))

        alignments = [alignment]
        readAlignments = ReadAlignments(Read('id0', 'aaa'), alignments)
        best = bestAlignment(readAlignments)
        self.assertEqual('Seq 1', best.hitTitle)
        self.assertEqual(44, best.hitLength)

    def testThreeAlignments(self):
        """
        When three alignments are present, the one with the lowest first HSP
        must be returned by bestAlignment.
        """
        alignment1 = Alignment(33, 'Seq 1')
        alignment1.addHsp(EValueHSP(score=10))
        alignment1.addHsp(EValueHSP(score=9))

        alignment2 = Alignment(44, 'Seq 2')
        alignment2.addHsp(EValueHSP(score=3))
        alignment2.addHsp(EValueHSP(score=2))

        alignment3 = Alignment(55, 'Seq 3')
        alignment3.addHsp(EValueHSP(score=20))
        alignment3.addHsp(EValueHSP(score=19))

        alignments = [alignment1, alignment2, alignment3]
        readAlignments = ReadAlignments(Read('id0', 'aaa'), alignments)
        best = bestAlignment(readAlignments)
        self.assertEqual('Seq 2', best.hitTitle)
        self.assertEqual(44, best.hitLength)

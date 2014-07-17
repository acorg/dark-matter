from unittest import TestCase
from json import dumps
from mock import patch

from dark import nicola
from dark.blast import BlastRecords
from mocking import mockOpen


class TestUtilitiesFunctions(TestCase):
    """
    Tests for the general utilities functions in nicola
    (_getBird, computePercentId, getPositionInMatrix).
    """

    def test_getBirdNoBird(self):
        title = "I'm a random title and don't contain a bird"
        result = nicola._getBird(title)
        self.assertEqual(result, None)

    def test_getBirdGull(self):
        title = "I'm a gull"
        result = nicola._getBird(title)
        self.assertEqual(result, 'gull')

    def test_getBirdDuck(self):
        title = "I'm a duck"
        result = nicola._getBird(title)
        self.assertEqual(result, 'duck')

    def testComputePercentIdIdentical(self):
        params = {
            'application': 'BLASTN',
        }
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 3.29804,
                            'sbjct': 'TACCCTGCGGCCCGCTACGGCTGG',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    "title": "Merkel1"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(params) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            records = list(BlastRecords('file.json').records(
                oneAlignmentPerRead=False))
            for alignment in records[0].alignments:
                result = nicola.computePercentId(alignment)
            self.assertEqual(result, 1.0)

    def testComputePercentIdDifferent(self):
        params = {
            'application': 'BLASTN',
        }
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 3.29804,
                            'sbjct': '------------------------',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    "title": "Merkel1"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(params) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            records = list(BlastRecords('file.json').records(
                oneAlignmentPerRead=False))
            for alignment in records[0].alignments:
                result = nicola.computePercentId(alignment)
            self.assertEqual(result, 0.0)

    def testComputePercentIdRealistic(self):
        params = {
            'application': 'BLASTN',
        }
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 3.29804,
                            'sbjct': 'TACCCTGCGGCCCGC-ACGGCTGG',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    "title": "Merkel1"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(params) + '\n' +
                              dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            records = list(BlastRecords('file.json').records(
                oneAlignmentPerRead=False))
            for alignment in records[0].alignments:
                result = nicola.computePercentId(alignment)
            self.assertEqual(result, 0.9583333333333334)

    def testGetPositionInMatrix(self):
        title = 'hi'
        matrix = [[0, 0, 0], [0, 'hi', 0], [0, 0, 0]]
        x, y = nicola.getPositionInMatrix(matrix, title)
        self.assertEqual([x, y], [1, 1])

    def testGetPositionInMatrixNotThere(self):
        title = 'hi'
        matrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        x, y = nicola.getPositionInMatrix(matrix, title)
        self.assertEqual([x, y], [False, False])

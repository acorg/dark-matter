from unittest import TestCase
from dark.smith_waterman import smith_waterman


class TestSmithWaterman(TestCase):
    """
    Test the smith_waterman function. 
    With match +1, mismatch -1 and gap -1.
    """

    def testZeroInput(self):
        seq1 = ''
        seq2 = ''
        result = smith_waterman(seq1, seq2)
        self.assertEqual(result, [])

    def testOneInput(self):
        seq1 = 'agtcagtcagtc'
        seq2 = ''
        result = smith_waterman(seq1, seq2)
        self.assertEqual(result, [])

    def testTwoInput(self):
        seq1 = 'agtcagtcagtc'
        seq2 = 'cgaatcg'
        result = smith_waterman(seq1, seq2)
        p = ['agtcagtcagtc', '      || |  ', '--cgaatc-g--']
        self.assertEqual(result, p)

    def testWikiAnswer(self):
        seq1 = 'ACACACTA'
        seq2 = 'AGCACACA'
        result = smith_waterman(seq1, seq2)
        p = ['-ACACACTA', '  ||||| |', 'AGCACAC-A']
        self.assertEqual(result, p)

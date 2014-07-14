from unittest import TestCase
from dark.smith_waterman import smith_waterman


class TestSmithWaterman(TestCase):
    """
    Test the smith_waterman function.
    With match +1, mismatch -1 and gap -1.
    """

    def testNonDNAString(self):
        """
        If the mismatch value passed is positive, an exception
        must be raised.
        """
        self.assertRaises(ValueError, smith_waterman, 'xxx', 'yyy')

    def testPositiveMismatch(self):
        """
        If the mismatch value passed is positive, an exception
        must be raised.
        """
        self.assertRaises(ValueError, smith_waterman, '', '', mismatch=3)

    def testZeroMismatch(self):
        """
        If the mismatch value passed is zero, an exception
        must be raised.
        """
        self.assertRaises(ValueError, smith_waterman, '', '', mismatch=0)

    def testPositiveGap(self):
        """
        If the gap value passed is positive, an exception
        must be raised.
        """
        self.assertRaises(ValueError, smith_waterman, '', '', gap=3)

    def testZeroGap(self):
        """
        If the gap value passed is zero, an exception
        must be raised.
        """
        self.assertRaises(ValueError, smith_waterman, '', '', gap=0)

    def testFirstSequenceEmpty(self):
        """
        If the first sequence passed is empty, the result should be a
        list of three empty strings.
        """
        seq1 = ''
        seq2 = 'agtcagtcagtc'
        result = smith_waterman(seq1, seq2)
        self.assertEqual(result, ['', '', ''])

    def testSecondSequenceEmpty(self):
        """
        If the second sequence passed is empty, the result should be a
        list of three empty strings.
        """
        seq1 = 'agtcagtcagtc'
        seq2 = ''
        result = smith_waterman(seq1, seq2)
        self.assertEqual(result, ['', '', ''])

    def testBothSequencesEmpty(self):
        """
        If two empty sequences are passed, the result should be a list of
        three empty strings.
        """
        seq1 = ''
        seq2 = ''
        result = smith_waterman(seq1, seq2)
        self.assertEqual(result, ['', '', ''])

    def testGapAtStartOfSeq1(self):
        seq1 = 'gaatcg'
        seq2 = 'cgaatcg'
        result = smith_waterman(seq1, seq2)
        p = ['-GAATCG',
             ' ||||||',
             'CGAATCG']
        self.assertEqual(result, p)

    def testGapAtStartOfSeq2(self):
        seq1 = 'cgaatcg'
        seq2 = 'gaatcg'
        result = smith_waterman(seq1, seq2)
        p = ['CGAATCG',
             ' ||||||',
             '-GAATCG']
        self.assertEqual(result, p)

    def testGapAtEndOfSeq1(self):
        seq1 = 'cgaatc'
        seq2 = 'cgaatcg'
        result = smith_waterman(seq1, seq2)
        p = ['CGAATC-',
             '|||||| ',
             'CGAATCG']
        self.assertEqual(result, p)

    def testGapAtEndOfSeq2(self):
        seq1 = 'cgaatcg'
        seq2 = 'cgaatc'
        result = smith_waterman(seq1, seq2)
        p = ['CGAATCG',
             '|||||| ',
             'CGAATC-']
        self.assertEqual(result, p)

    def testGapAtBothEndsOfSeq1(self):
        seq1 = 'gaatc'
        seq2 = 'cgaatcg'
        result = smith_waterman(seq1, seq2)
        p = ['-GAATC-',
             ' ||||| ',
             'CGAATCG']
        self.assertEqual(result, p)

    def testGapAtBothEndsOfSeq2(self):
        seq1 = 'cgaatcg'
        seq2 = 'gaatc'
        result = smith_waterman(seq1, seq2)
        p = ['CGAATCG',
             ' ||||| ',
             '-GAATC-']
        self.assertEqual(result, p)

    def testAlignmentWithGapInMiddle(self):
        seq1 = 'agtcagtcagtc'
        seq2 = 'cgaatcg'
        result = smith_waterman(seq1, seq2)
        p = ['AGTCAGTCAGTC', '      || |  ', '--CGAATC-G--']
        self.assertEqual(result, p)

    def testTwoEqualSequences(self):
        seq1 = seq2 = 'cgaatcg'
        result = smith_waterman(seq1, seq2)
        p = ['CGAATCG',
             '|||||||',
             'CGAATCG']
        self.assertEqual(result, p)

    def testTwoCompletelyDifferentSequences(self):
        """
        When two completely different sequences are given, the result
        should be the two sequences with an empty alignment.
        """
        seq1 = 'aaaaaa'
        seq2 = 'gggggg'
        result = smith_waterman(seq1, seq2)
        p = ['AAAAAA',
             '      ',
             'GGGGGG']
        self.assertEqual(result, p)

    def testWikiAnswer(self):
        """
        Test the example given in Wikipedia:
        http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        """
        seq1 = 'ACACACTA'
        seq2 = 'AGCACACA'
        result = smith_waterman(seq1, seq2, match=2)
        p = ['A-CACACTA',
             '| ||||| |',
             'AGCACAC-A']
        self.assertEqual(result, p)

    def testWikiAnswerWithMatchOne(self):
        """
        Test the example given in Wikipedia
        http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        Wikipedia uses a match score of two, here we use a score of one.
        """
        seq1 = 'ACACACTA'
        seq2 = 'AGCACACA'
        result = smith_waterman(seq1, seq2, match=1)
        p = ['-ACACACTA', '  ||||| |', 'AGCACAC-A']
        self.assertEqual(result, p)

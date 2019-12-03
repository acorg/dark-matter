from unittest import TestCase

from dark.hsp import HSP, LSP


class TestHSP(TestCase):
    """
    Tests of the L{dark.hsp.HSP} class.
    """

    def testExpectedAttributes(self):
        """
        An HSP must have the expected attributes.
        """
        hsp = HSP(7, readStart=1, readEnd=2,
                  readStartInSubject=3, readEndInSubject=4,
                  subjectStart=5, subjectEnd=6,
                  readMatchedSequence='aaa', subjectMatchedSequence='ccc',
                  readFrame=8, subjectFrame=9)
        self.assertEqual(1, hsp.readStart)
        self.assertEqual(2, hsp.readEnd)
        self.assertEqual(3, hsp.readStartInSubject)
        self.assertEqual(4, hsp.readEndInSubject)
        self.assertEqual(5, hsp.subjectStart)
        self.assertEqual(6, hsp.subjectEnd)
        self.assertEqual('aaa', hsp.readMatchedSequence)
        self.assertEqual('ccc', hsp.subjectMatchedSequence)
        self.assertEqual(7, hsp.score.score)
        self.assertEqual(8, hsp.readFrame)
        self.assertEqual(9, hsp.subjectFrame)

    def testEqual(self):
        """
        Two HSPs must compare properly with ==
        """
        self.assertEqual(HSP(7), HSP(7))

    def testLt(self):
        """
        Two HSPs must compare properly with <
        """
        self.assertTrue(HSP(7) < HSP(8))

    def testBetterThanTrue(self):
        """
        Two HSP scores must compare properly with betterThan if the passed
        score is worse than the score of the HSP.
        """
        self.assertTrue(HSP(7).betterThan(5))

    def testBetterThanFalse(self):
        """
        Two HSP scores must compare properly with betterThan if the passed
        score is better than the score of the HSP.
        """
        self.assertFalse(HSP(5).betterThan(7))

    def testToDict(self):
        """
        The toDict method must return the expected dictionary.
        """
        hsp = HSP(0, readStart=1, readEnd=2,
                  readStartInSubject=3, readEndInSubject=4,
                  subjectStart=5, subjectEnd=6,
                  readMatchedSequence='aaa', subjectMatchedSequence='ccc',
                  readFrame=7, subjectFrame=8, identicalCount=9,
                  percentIdentical=99.3, positiveCount=10, percentPositive=3.0)

        self.assertEqual(
            {
                'score': 0,
                'readStart': 1,
                'readEnd': 2,
                'readStartInSubject': 3,
                'readEndInSubject': 4,
                'subjectStart': 5,
                'subjectEnd': 6,
                'readFrame': 7,
                'subjectFrame': 8,
                'identicalCount': 9,
                'percentIdentical': 99.3,
                'positiveCount': 10,
                'percentPositive': 3.0,
                'readMatchedSequence': 'aaa',
                'subjectMatchedSequence': 'ccc',
            },
            hsp.toDict())


class TestLSP(TestCase):
    """
    Tests of the L{dark.hsp.blast.LSP} class.
    """

    def testExpectedAttributes(self):
        """
        An LSP must have the expected attributes.
        """
        lsp = LSP(7, readStart=1, readEnd=2,
                  readStartInSubject=3, readEndInSubject=4,
                  subjectStart=5, subjectEnd=6,
                  readMatchedSequence='aaa', subjectMatchedSequence='ccc',
                  readFrame=8, subjectFrame=9)
        self.assertEqual(1, lsp.readStart)
        self.assertEqual(2, lsp.readEnd)
        self.assertEqual(3, lsp.readStartInSubject)
        self.assertEqual(4, lsp.readEndInSubject)
        self.assertEqual(5, lsp.subjectStart)
        self.assertEqual(6, lsp.subjectEnd)
        self.assertEqual('aaa', lsp.readMatchedSequence)
        self.assertEqual('ccc', lsp.subjectMatchedSequence)
        self.assertEqual(7, lsp.score.score)
        self.assertEqual(8, lsp.readFrame)
        self.assertEqual(9, lsp.subjectFrame)

    def testEqual(self):
        """
        Two LSPs must compare properly with ==
        """
        self.assertEqual(LSP(7), LSP(7))

    def testLt(self):
        """
        Two LSPs must compare properly with <
        """
        self.assertTrue(LSP(8) < LSP(7))

    def testBetterThanTrue(self):
        """
        Two LSP scores must compare properly with betterThan if the passed
        score is worse than the score of the LSP.
        """
        self.assertTrue(LSP(5).betterThan(7))

    def testBetterThanFalse(self):
        """
        Two LSP scores must compare properly with betterThan if the passed
        score is better than the score of the LSP.
        """
        self.assertFalse(LSP(7).betterThan(5))

    def testToDict(self):
        """
        The toDict method must return the expected dictionary.
        """
        lsp = LSP(0, readStart=1, readEnd=2,
                  readStartInSubject=3, readEndInSubject=4,
                  subjectStart=5, subjectEnd=6,
                  readMatchedSequence='aaa', subjectMatchedSequence='ccc',
                  readFrame=7, subjectFrame=8, identicalCount=9,
                  percentIdentical=99.3, positiveCount=10,
                  percentPositive=9.9)

        self.assertEqual(
            {
                'score': 0,
                'readStart': 1,
                'readEnd': 2,
                'readStartInSubject': 3,
                'readEndInSubject': 4,
                'subjectStart': 5,
                'subjectEnd': 6,
                'readFrame': 7,
                'subjectFrame': 8,
                'identicalCount': 9,
                'percentIdentical': 99.3,
                'positiveCount': 10,
                'percentPositive': 9.9,
                'readMatchedSequence': 'aaa',
                'subjectMatchedSequence': 'ccc',
            },
            lsp.toDict())

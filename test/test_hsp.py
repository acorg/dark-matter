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
                  readMatchedSequence='aaa', subjectMatchedSequence='ccc')
        self.assertEqual(1, hsp.readStart)
        self.assertEqual(2, hsp.readEnd)
        self.assertEqual(3, hsp.readStartInSubject)
        self.assertEqual(4, hsp.readEndInSubject)
        self.assertEqual(5, hsp.subjectStart)
        self.assertEqual(6, hsp.subjectEnd)
        self.assertEqual('aaa', hsp.readMatchedSequence)
        self.assertEqual('ccc', hsp.subjectMatchedSequence)
        self.assertEqual(7, hsp.score.score)

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


class TestLSP(TestCase):
    """
    Tests of the L{dark.hsp.blast.LSP} class.
    """

    def testExpectedAttributes(self):
        """
        An LSP must have the expected attributes.
        """
        hsp = LSP(7, readStart=1, readEnd=2,
                  readStartInSubject=3, readEndInSubject=4,
                  subjectStart=5, subjectEnd=6,
                  readMatchedSequence='aaa', subjectMatchedSequence='ccc')
        self.assertEqual(1, hsp.readStart)
        self.assertEqual(2, hsp.readEnd)
        self.assertEqual(3, hsp.readStartInSubject)
        self.assertEqual(4, hsp.readEndInSubject)
        self.assertEqual(5, hsp.subjectStart)
        self.assertEqual(6, hsp.subjectEnd)
        self.assertEqual('aaa', hsp.readMatchedSequence)
        self.assertEqual('ccc', hsp.subjectMatchedSequence)
        self.assertEqual(7, hsp.score.score)

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

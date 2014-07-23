from unittest import TestCase

from dark.hsp import HSP


class TestHSP(TestCase):
    """
    Tests of the L{dark.hsp.HSP} class.
    """

    def testExpectedAttributes(self):
        """
        An HSP must have the expected attributes.
        """
        hsp = HSP(readStart=1, readEnd=2,
                  readStartInHit=3, readEndInHit=4,
                  hitStart=5, hitEnd=6,
                  readMatchedSequence='aaa', hitMatchedSequence='ccc',
                  score=7)
        self.assertEqual(1, hsp.readStart)
        self.assertEqual(2, hsp.readEnd)
        self.assertEqual(3, hsp.readStartInHit)
        self.assertEqual(4, hsp.readEndInHit)
        self.assertEqual(5, hsp.hitStart)
        self.assertEqual(6, hsp.hitEnd)
        self.assertEqual('aaa', hsp.readMatchedSequence)
        self.assertEqual('ccc', hsp.hitMatchedSequence)
        self.assertEqual(7, hsp.score)

    def testEq(self):
        """
        Two HSPs must compare properly with ==
        """
        hsp1 = HSP(score=7)
        hsp2 = HSP(score=7)
        self.assertEqual(hsp1, hsp2)

    def testLt(self):
        """
        Two HSPs must compare properly with <
        """
        hsp1 = HSP(score=7)
        hsp2 = HSP(score=8)
        self.assertTrue(hsp1 < hsp2)

    def testBetterThanTrue(self):
        """
        Two HSP scores must compare properly with betterThan if the passed
        score is worse than the score of the HSP.
        """
        self.assertTrue(HSP(score=7).betterThan(5))

    def testBetterThanFalse(self):
        """
        Two HSP scores must compare properly with betterThan if the passed
        score is better than the score of the HSP.
        """
        self.assertFalse(HSP(score=5).betterThan(7))

from unittest import TestCase

from dark.gor4 import GOR4


class TestGOR4(TestCase):
    """
    Tests for the L{dark.gor4.GOR4} class.
    """

    def testNothing(self):
        """
        """
        gor4 = GOR4()
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        expected = 'CCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEECEEC'
        result = gor4.predict(seq)
        self.assertEqual(expected, result['predictions'])

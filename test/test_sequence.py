from unittest import TestCase
from dark.sequence import findPrimer, findPrimerBidi, findPrimerBidiLimits
from Bio.Seq import Seq


class TestFindPrimer(TestCase):
    """
    Tests for the dark.sequence.findPrimer function.
    """

    def testNotFound(self):
        """
        If a primer is not found, the empty list must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual([], findPrimer("BLAH", seq))

    def testFoundAtStart(self):
        """
        If a primer is found at the start of a sequence, a list containing 0
        must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual([0], findPrimer("AC", seq))

    def testFoundAtEnd(self):
        """
        If a primer is found at the end of a sequence, the correct value
        must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual([2], findPrimer("GT", seq))

    def testFoundMultiple(self):
        """
        If a primer is found multiple times, the correct value
        must be returned.
        """
        seq = Seq("ACGTACGT")
        self.assertEqual([0, 4], findPrimer("ACG", seq))

    def testOverlapping(self):
        """
        If a primer is present twice but is overlapping, only the first
        instance should be returned.
        """
        seq = Seq("GAAA")
        self.assertEqual([1], findPrimer("AA", seq))


class TestFindPrimerBidi(TestCase):
    """
    Tests for the dark.sequence.findPrimerBidi function.
    """

    def testNotFound(self):
        """
        If a primer is not found, empty lists must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual(([], []), findPrimerBidi("BLAH", seq))

    def testFoundStartEnd(self):
        """
        If a primer is found in both directions in a sequence (start of
        the forward sequence, end of the reverse complement), the
        correct value must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual(([0], [2]), findPrimerBidi("AC", seq))

    def testFoundEndStart(self):
        """
        If a primer is found in both directions in a sequence (end of
        the forward sequence, start of the reverse complement), the
        correct value must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual(([2], [0]), findPrimerBidi("GT", seq))

    def testFoundMultiple(self):
        """
        If a primer is found multiple times, the correct value
        must be returned.
        """
        seq = Seq("ACGTACGT")
        self.assertEqual(([0, 4], [1, 5]), findPrimerBidi("ACG", seq))

    def testOverlappingForwards(self):
        """
        If a primer is present twice forwards but is overlapping, only
        the first instance should be returned.
        """
        seq = Seq("GAAA")
        self.assertEqual(([1], []), findPrimerBidi("AA", seq))

    def testOverlappingBackwards(self):
        """
        If a primer is present twice backwards but is overlapping, only
        the first instance should be returned.
        """
        seq = Seq("GTTT")
        self.assertEqual(([], [1]), findPrimerBidi("AA", seq))


class TestFindPrimerBidiLimits(TestCase):
    """
    Tests for the dark.sequence.findPrimerBidiLimits function.
    """

    def testNotFound(self):
        """
        If a primer is not found, the returned offsets must include
        the whole sequence.
        """
        seq = Seq("ACGT")
        self.assertEqual((0, 4), findPrimerBidiLimits("BLAH", seq))

    def testFoundStartEnd(self):
        """
        If a primer is found in both directions in a sequence (start of
        the forward sequence, end of the reverse complement), the
        correct value must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual((2, 2), findPrimerBidiLimits("AC", seq))

    def testFoundEndStart(self):
        """
        If a primer is found in both directions in a sequence (end of
        the forward sequence, start of the reverse complement), the
        correct value must be returned.
        """
        seq = Seq("ACGT")
        self.assertEqual((4, 4), findPrimerBidiLimits("GT", seq))

    def testFoundMultiple(self):
        """
        If a primer is found multiple times, the correct value
        must be returned.
        """
        seq = Seq("ACGTACGT")
        self.assertEqual((7, 8), findPrimerBidiLimits("ACG", seq))

    def testOverlappingForwards(self):
        """
        If a primer is present twice forwards but is overlapping, only
        the first instance should be returned.
        """
        seq = Seq("GAAA")
        self.assertEqual((3, 4), findPrimerBidiLimits("AA", seq))

    def testOverlappingBackwards(self):
        """
        If a primer is present twice backwards but is overlapping, only
        the first instance should be returned.
        """
        seq = Seq("GTTT")
        self.assertEqual((0, 1), findPrimerBidiLimits("AA", seq))

    def testLonger(self):
        """
        Test a longer sequence.
        """
        seq = Seq("AAAAAAAAAA" "GGGGGGGGGG" "AAAAAAAAAA" "AAAAAAAAAA")
        self.assertEqual((20, 40), findPrimerBidiLimits("GGGGGGGGGG", seq))

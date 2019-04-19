from unittest import TestCase
from six import assertRaisesRegex

from dark.cigar import CINS, CDEL, CMATCH, CEQUAL, CDIFF, dna2cigar


class TestConstants(TestCase):
    """
    Test constants in dark.cigar
    """
    def testOperations(self):
        """
        Make sure the CIGAR operations have the expected one-letter codes.
        """
        self.assertEqual('I', CINS)
        self.assertEqual('D', CDEL)
        self.assertEqual('M', CMATCH)
        self.assertEqual('=', CEQUAL)
        self.assertEqual('X', CDIFF)


class TestDNA2CIGAR(TestCase):
    """
    Test the dna2cigar function.
    """
    def testUnequalStrings(self):
        """
        If two unequal strings are passed, a ValueError must be raised.
        """
        error = r"^Sequences 'hey' and 'there' of unequal length \(3 != 5\)\.$"
        assertRaisesRegex(self, ValueError, error, dna2cigar, 'hey', 'there')

    def testEmptyStrings(self):
        """
        If two empty strings are passed, a ValueError must be raised.
        """
        error = r"^Two sequences of zero length were passed\.$"
        assertRaisesRegex(self, ValueError, error, dna2cigar, '', '')

    def testEqualStringsLength1(self):
        """
        If two equal strings of length 1 are passed, the correct CIGAR string
        must be returned.
        """
        self.assertEqual('1=', dna2cigar('A', 'A'))

    def testEqualStringsLength1MismatchedCase(self):
        """
        If two strings of length 1 are passed and the're the same but for
        case, the CIGAR string must indicate that they do not match.
        """
        self.assertEqual('1X', dna2cigar('A', 'a'))

    def testEqualStringsLength1Concise(self):
        """
        If two equal strings of length 1 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual('1M', dna2cigar('A', 'A', concise=True))

    def testUnequalStringsLength1(self):
        """
        If two unequal strings of length 1 are passed, the correct CIGAR
        string must be returned.
        """
        self.assertEqual('1X', dna2cigar('A', 'G'))

    def testUnequalStringsLength1Concise(self):
        """
        If two unequal strings of length 1 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual('1M', dna2cigar('A', 'G', concise=True))

    def testEqualStringsLength2(self):
        """
        If two equal strings of length 2 are passed, the correct CIGAR string
        must be returned.
        """
        self.assertEqual('2=', dna2cigar('AT', 'AT'))

    def testEqualStringsLength2Concise(self):
        """
        If two equal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual('2M', dna2cigar('AT', 'AT', concise=True))

    def testUnequalStringsLength2(self):
        """
        If two unequal strings of length 2 are passed, the correct CIGAR
        string must be returned.
        """
        self.assertEqual('2X', dna2cigar('AT', 'GC'))

    def testUnequalStringsLength2Concise(self):
        """
        If two unequal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual('2M', dna2cigar('AT', 'GC', concise=True))

    def testMixedEqualityStringsLength2(self):
        """
        If two strings of length 2 are passed and they match in just one
        place, the correct CIGAR string must be returned.
        """
        self.assertEqual('1=1X', dna2cigar('AG', 'AT'))

    def testMixedEqualityStringsLength2Concise(self):
        """
        If two equal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual('2M', dna2cigar('AG', 'AT', concise=True))

    def testMixedEqualityStringsLength2Alt(self):
        """
        If two strings of length 2 are passed and they match in just one
        place, the correct CIGAR string must be returned.
        """
        self.assertEqual('1X1=', dna2cigar('GA', 'TA'))

    def testMixedEqualityStringsLength2ConciseAlt(self):
        """
        If two equal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual('2M', dna2cigar('GA', 'TA', concise=True))

    def testMixedMatch(self):
        """
        If two strings with matching and non-matching regions must result
        in the expected CIGAR string.
        """
        self.assertEqual('3=4X5=2X',
                         dna2cigar('ACGTTTTCCCTTGG',
                                   'ACGAAAACCCTTTT'))

    def testMixedMatchConcise(self):
        """
        If two strings with matching and non-matching regions must result
        in the expected CIGAR string when the concise argument is True.
        """
        self.assertEqual('14M',
                         dna2cigar('ACGTTTTCCCTTGG',
                                   'ACGAAAACCCTTTT', concise=True))

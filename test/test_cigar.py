from unittest import TestCase

from dark.cigar import CINS, CDEL, CMATCH, CEQUAL, CDIFF, dna2cigar, makeCigar


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
        self.assertRaisesRegex(ValueError, error, dna2cigar, 'hey', 'there')

    def testEmptyStrings(self):
        """
        If two empty strings are passed, a ValueError must be raised.
        """
        error = r"^Two sequences of zero length were passed\.$"
        self.assertRaisesRegex(ValueError, error, dna2cigar, '', '')

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


class TestMakeCigar(TestCase):
    """
    Test making CIGAR strings.
    """
    def testEmptyReference(self):
        """
        If the reference is empty, a ValueError must be raised.
        """
        error = r"^Empty reference$"
        self.assertRaisesRegex(ValueError, error, makeCigar, '', 'ACGT')

    def testEmptyQuery(self):
        """
        If the query is empty, a ValueError must be raised.
        """
        error = r"^Empty query$"
        self.assertRaisesRegex(ValueError, error, makeCigar, 'ACGT', '')

    def testOneBaseInitialMatch(self):
        """
        If the query has one base that matches the start of the reference, 1M
        should be the result.
        """
        self.assertEqual('1M', makeCigar('ACGT',
                                         'A'))

    def testTwoBasesInitialMatch(self):
        """
        If the query has two bases that matches the start of the reference,
        1M1M should be the result.
        """
        self.assertEqual('2M', makeCigar('ACGT',
                                         'AT'))

    def testOneBaseFinalMatch(self):
        """
        If the query has one base that matches the end of the reference, 1M
        should be the result.
        """
        self.assertEqual('1M', makeCigar('ACGT',
                                         '   A'))

    def testTwoBasesFinalMatch(self):
        """
        If the query has two bases that matches the end of the reference, 1M1M
        should be the result.
        """
        self.assertEqual('2M', makeCigar('ACGT',
                                         '  AT'))

    def testEqualStrings(self):
        """
        If the query exactly overlaps the reference we should get 'M' as many
        times as the length of the sequences.
        """
        self.assertEqual('4M', makeCigar('ACGT',
                                         'ATGG'))

    def testLeftClippingOne(self):
        """
        Test that left soft clipping of one base works.
        """
        self.assertEqual('2S2M', makeCigar('  ACGT',
                                           'ATAT'))

    def testLeftClippingTwo(self):
        """
        Test that left soft clipping of two bases works.
        """
        self.assertEqual('2S1M', makeCigar('  ACGT',
                                           'TCA'))

    def testRightClippingUnpaddedOne(self):
        """
        Test that right soft clipping of one base works when the reference
        does not have enough padding to match the length of the query.
        """
        self.assertEqual('2M1S', makeCigar('  ACGT',
                                           '    ATA'))

    def testRightClippingUnpaddedTwo(self):
        """
        Test that right soft clipping of two bases works when the reference
        does not have enough padding to match the length of the query.
        """
        self.assertEqual('2M2S', makeCigar('  ACGT',
                                           '    ATAA'))

    def testRightClippingPaddedOne(self):
        """
        Test that right soft clipping of one base works when the reference
        is padded on the right.
        """
        self.assertEqual('2M1S', makeCigar('  ACGT ',
                                           '    ATA'))

    def testRightClippingPaddedTwo(self):
        """
        Test that right soft clipping of two bases works when the reference
        is padded on the right.
        """
        self.assertEqual('2M2S', makeCigar('  ACGT  ',
                                           '    ATAA'))

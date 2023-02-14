from unittest import TestCase
from pysam import CSOFT_CLIP, CHARD_CLIP, CMATCH, CINS

from dark.cigar import (
    CINS_STR,
    CDEL_STR,
    CMATCH_STR,
    CEQUAL_STR,
    CDIFF_STR,
    dna2cigar,
    makeCigar,
    cigarTuplesToOperations,
    softClippedOffset,
    insertionOffset,
)


class TestConstants(TestCase):
    """
    Test constants in dark.cigar
    """

    def testOperations(self):
        """
        Make sure the CIGAR operation strings have the expected one-letter
        codes.
        """
        self.assertEqual("I", CINS_STR)
        self.assertEqual("D", CDEL_STR)
        self.assertEqual("M", CMATCH_STR)
        self.assertEqual("=", CEQUAL_STR)
        self.assertEqual("X", CDIFF_STR)


class TestDNA2CIGAR(TestCase):
    """
    Test the dna2cigar function.
    """

    def testUnequalStrings(self):
        """
        If two unequal strings are passed, a ValueError must be raised.
        """
        error = r"^Sequences 'hey' and 'there' of unequal length \(3 != 5\)\.$"
        self.assertRaisesRegex(ValueError, error, dna2cigar, "hey", "there")

    def testEmptyStrings(self):
        """
        If two empty strings are passed, a ValueError must be raised.
        """
        error = r"^Two sequences of zero length were passed\.$"
        self.assertRaisesRegex(ValueError, error, dna2cigar, "", "")

    def testEqualStringsLength1(self):
        """
        If two equal strings of length 1 are passed, the correct CIGAR string
        must be returned.
        """
        self.assertEqual("1=", dna2cigar("A", "A"))

    def testEqualStringsLength1MismatchedCase(self):
        """
        If two strings of length 1 are passed and the're the same but for
        case, the CIGAR string must indicate that they do not match.
        """
        self.assertEqual("1X", dna2cigar("A", "a"))

    def testEqualStringsLength1Concise(self):
        """
        If two equal strings of length 1 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual("1M", dna2cigar("A", "A", concise=True))

    def testUnequalStringsLength1(self):
        """
        If two unequal strings of length 1 are passed, the correct CIGAR
        string must be returned.
        """
        self.assertEqual("1X", dna2cigar("A", "G"))

    def testUnequalStringsLength1Concise(self):
        """
        If two unequal strings of length 1 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual("1M", dna2cigar("A", "G", concise=True))

    def testEqualStringsLength2(self):
        """
        If two equal strings of length 2 are passed, the correct CIGAR string
        must be returned.
        """
        self.assertEqual("2=", dna2cigar("AT", "AT"))

    def testEqualStringsLength2Concise(self):
        """
        If two equal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual("2M", dna2cigar("AT", "AT", concise=True))

    def testUnequalStringsLength2(self):
        """
        If two unequal strings of length 2 are passed, the correct CIGAR
        string must be returned.
        """
        self.assertEqual("2X", dna2cigar("AT", "GC"))

    def testUnequalStringsLength2Concise(self):
        """
        If two unequal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual("2M", dna2cigar("AT", "GC", concise=True))

    def testMixedEqualityStringsLength2(self):
        """
        If two strings of length 2 are passed and they match in just one
        place, the correct CIGAR string must be returned.
        """
        self.assertEqual("1=1X", dna2cigar("AG", "AT"))

    def testMixedEqualityStringsLength2Concise(self):
        """
        If two equal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual("2M", dna2cigar("AG", "AT", concise=True))

    def testMixedEqualityStringsLength2Alt(self):
        """
        If two strings of length 2 are passed and they match in just one
        place, the correct CIGAR string must be returned.
        """
        self.assertEqual("1X1=", dna2cigar("GA", "TA"))

    def testMixedEqualityStringsLength2ConciseAlt(self):
        """
        If two equal strings of length 2 are passed and the concise argument
        is True, the correct CIGAR string must be returned.
        """
        self.assertEqual("2M", dna2cigar("GA", "TA", concise=True))

    def testMixedMatch(self):
        """
        If two strings with matching and non-matching regions must result
        in the expected CIGAR string.
        """
        self.assertEqual("3=4X5=2X", dna2cigar("ACGTTTTCCCTTGG", "ACGAAAACCCTTTT"))

    def testMixedMatchConcise(self):
        """
        If two strings with matching and non-matching regions must result
        in the expected CIGAR string when the concise argument is True.
        """
        self.assertEqual(
            "14M", dna2cigar("ACGTTTTCCCTTGG", "ACGAAAACCCTTTT", concise=True)
        )


class TestMakeCigar(TestCase):
    """
    Test making CIGAR strings.
    """

    def testEmptyReference(self):
        """
        If the reference is empty, a ValueError must be raised.
        """
        error = r"^Empty reference$"
        self.assertRaisesRegex(ValueError, error, makeCigar, "", "ACGT")

    def testEmptyQuery(self):
        """
        If the query is empty, a ValueError must be raised.
        """
        error = r"^Empty query$"
        self.assertRaisesRegex(ValueError, error, makeCigar, "ACGT", "")

    def testOneBaseInitialMatch(self):
        """
        If the query has one base that matches the start of the reference, 1M
        should be the result.
        """
        self.assertEqual("1M", makeCigar("ACGT", "A"))

    def testTwoBasesInitialMatch(self):
        """
        If the query has two bases that matches the start of the reference,
        1M1M should be the result.
        """
        self.assertEqual("2M", makeCigar("ACGT", "AT"))

    def testOneBaseFinalMatch(self):
        """
        If the query has one base that matches the end of the reference, 1M
        should be the result.
        """
        self.assertEqual("1M", makeCigar("ACGT", "   A"))

    def testTwoBasesFinalMatch(self):
        """
        If the query has two bases that matches the end of the reference, 1M1M
        should be the result.
        """
        self.assertEqual("2M", makeCigar("ACGT", "  AT"))

    def testEqualStrings(self):
        """
        If the query exactly overlaps the reference we should get 'M' as many
        times as the length of the sequences.
        """
        self.assertEqual("4M", makeCigar("ACGT", "ATGG"))

    def testLeftClippingOne(self):
        """
        Test that left soft clipping of one base works.
        """
        self.assertEqual("1S2M", makeCigar(" ACGT", "TAT"))

    def testLeftClippingTwo(self):
        """
        Test that left soft clipping of two bases works.
        """
        self.assertEqual("2S1M", makeCigar("  ACGT", "TCA"))

    def testRightClippingUnpaddedOne(self):
        """
        Test that right soft clipping of one base works when the reference
        does not have enough padding to match the length of the query.
        """
        self.assertEqual("2M1S", makeCigar("  ACGT", "    ATA"))

    def testRightClippingUnpaddedTwo(self):
        """
        Test that right soft clipping of two bases works when the reference
        does not have enough padding to match the length of the query.
        """
        self.assertEqual("2M2S", makeCigar("  ACGT", "    ATAA"))

    def testRightClippingPaddedOne(self):
        """
        Test that right soft clipping of one base works when the reference
        is padded on the right.
        """
        self.assertEqual("2M1S", makeCigar("  ACGT ", "    ATA"))

    def testRightClippingPaddedTwo(self):
        """
        Test that right soft clipping of two bases works when the reference
        is padded on the right.
        """
        self.assertEqual("2M2S", makeCigar("  ACGT  ", "    ATAA"))

    def testLeftInsert(self):
        """
        Test that an insert on the left works.
        """
        self.assertEqual("2I2M", makeCigar("--ACGT  ", "ATAC", noEdgeInsertions=False))

        # When edge insertions are not allowed, the first two bases are
        # soft-clipped (as in tests above).
        self.assertEqual("2S2M", makeCigar("--ACGT  ", "ATAC"))


class TestCigarTuplesToOperations(TestCase):
    """
    Test the cigarTuplesToOperations function.
    """

    def testEmpty(self):
        """
        If there are no tuples, the empty list must be returned.
        """
        self.assertEqual([], list(cigarTuplesToOperations([])))

    def testOneOfOne(self):
        """
        Test one tuple with one operation of length one.
        """
        self.assertEqual([CINS], list(cigarTuplesToOperations([(CINS, 1)])))

    def testOneOfTwo(self):
        """
        Test one tuple with one operation of length two.
        """
        self.assertEqual([CINS, CINS], list(cigarTuplesToOperations([(CINS, 2)])))

    def testTwo(self):
        """
        Test two tuples.
        """
        self.assertEqual(
            [CINS, CINS, CMATCH, CMATCH, CMATCH],
            list(cigarTuplesToOperations([(CINS, 2), (CMATCH, 3)])),
        )

    def testHardClippingIncluded(self):
        """
        Hard clipping operations should be included.
        """
        self.assertEqual(
            [CHARD_CLIP, CHARD_CLIP, CINS, CINS, CMATCH, CMATCH, CMATCH],
            list(cigarTuplesToOperations([(CHARD_CLIP, 2), (CINS, 2), (CMATCH, 3)])),
        )

    def testSkipHardClipping(self):
        """
        It must be possible to have hard clipping operations skipped.
        """
        self.assertEqual(
            [CINS, CINS, CMATCH, CMATCH, CMATCH],
            list(
                cigarTuplesToOperations(
                    [(CHARD_CLIP, 2), (CINS, 2), (CMATCH, 3)], includeHardClip=False
                )
            ),
        )


class TestSoftClippedOffset(TestCase):
    """
    Test the softClippedOffset function.
    """

    def testNoNonClippedBases(self):
        """
        If there are no preceding or following non-clipped bases, a ValueError
        must be raised.
        """
        error = (
            r"^Soft-clipped base with no following or preceding "
            r"non-hard-clipped bases\.$"
        )
        self.assertRaisesRegex(
            ValueError, error, softClippedOffset, 0, ((0, None),), (CSOFT_CLIP,)
        )

    def testMatchOneBefore(self):
        """
        Test that a soft-clipped base one site before a non-soft-clipped site
        returns the correct offset.
        """
        self.assertEqual(
            9, softClippedOffset(0, ((0, None), (1, 10)), (CSOFT_CLIP, CMATCH))
        )

    def testMatchTwoBefore(self):
        """
        Test that a soft-clipped base two sites before a non-soft-clipped site
        returns the correct offset.
        """
        self.assertEqual(
            8,
            softClippedOffset(
                0, ((0, None), (1, None), (2, 10)), (CSOFT_CLIP, CSOFT_CLIP, CMATCH)
            ),
        )

    def testMatchTwoBeforeWithHardClips(self):
        """
        Test that a soft-clipped base two sites before a non-soft-clipped site
        returns the correct offset.
        """
        self.assertEqual(
            8,
            softClippedOffset(
                2,
                ((0, None), (1, None), (2, None), (3, None), (4, 10)),
                (CHARD_CLIP, CHARD_CLIP, CSOFT_CLIP, CSOFT_CLIP, CMATCH),
            ),
        )

    def testMatchOneAfter(self):
        """
        Test that a soft-clipped base one site after a non-soft-clipped site
        returns the correct offset.
        """
        self.assertEqual(
            11, softClippedOffset(1, ((0, 10), (1, None)), (CMATCH, CSOFT_CLIP))
        )

    def testMatchTwoAfter(self):
        """
        Test that a soft-clipped base two sites after a non-soft-clipped site
        returns the correct offset.
        """
        self.assertEqual(
            12,
            softClippedOffset(
                2, ((0, 10), (1, None), (2, None)), (CMATCH, CSOFT_CLIP, CSOFT_CLIP)
            ),
        )

    def testMatchTwoAfterThenHardClips(self):
        """
        Test that a soft-clipped base two sites after a non-soft-clipped site
        returns the correct offset, including when there are also hard clips.
        """
        self.assertEqual(
            12,
            softClippedOffset(
                2,
                ((0, 10), (1, None), (2, None), (3, None), (4, None)),
                (CMATCH, CSOFT_CLIP, CSOFT_CLIP, CHARD_CLIP, CHARD_CLIP),
            ),
        )


class TestInsertionOffset(TestCase):
    """
    Test the insertionOffset function.
    """

    def testNoNonInsertionBases(self):
        """
        If there are no preceding or following non-inserted bases, a ValueError
        must be raised.
        """
        error = r"^Inserted base with no following or preceding reference " r"bases\.$"
        self.assertRaisesRegex(
            ValueError, error, insertionOffset, 0, ((0, None),), (CINS,)
        )

    def testMatchOneBefore(self):
        """
        Test that an inserted base one site before a non-inserted site
        returns the correct offset.
        """
        self.assertEqual(
            (False, 10), insertionOffset(0, ((0, None), (1, 10)), (CINS, CMATCH))
        )

    def testMatchTwoBefore(self):
        """
        Test that an inserted base two sites before a non-inserted site
        returns the correct offset.
        """
        self.assertEqual(
            (False, 10),
            insertionOffset(0, ((0, None), (1, None), (2, 10)), (CINS, CINS, CMATCH)),
        )

    def testMatchTwoBeforeWithHardClips(self):
        """
        Test that an inserted base two sites before a non-inserted site
        returns the correct offset.
        """
        self.assertEqual(
            (False, 10),
            insertionOffset(
                2,
                ((0, None), (1, None), (2, None), (3, None), (4, 10)),
                (CHARD_CLIP, CHARD_CLIP, CINS, CINS, CMATCH),
            ),
        )

    def testMatchOneAfter(self):
        """
        Test that an inserted base one site after a non-inserted site
        returns the correct offset.
        """
        self.assertEqual(
            (True, 11), insertionOffset(1, ((0, 10), (1, None)), (CMATCH, CINS))
        )

    def testMatchTwoAfter(self):
        """
        Test that an inserted base two sites after a non-inserted site
        returns the correct offset.
        """
        self.assertEqual(
            (True, 11),
            insertionOffset(2, ((0, 10), (1, None), (2, None)), (CMATCH, CINS, CINS)),
        )

    def testMatchTwoAfterThenHardClips(self):
        """
        Test that an inserted base two sites after a non-inserted site
        returns the correct offset, including when there are also hard clips.
        """
        self.assertEqual(
            (True, 11),
            insertionOffset(
                2,
                ((0, 10), (1, None), (2, None), (3, None), (4, None)),
                (CMATCH, CINS, CINS, CHARD_CLIP, CHARD_CLIP),
            ),
        )

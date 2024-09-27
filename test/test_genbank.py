from unittest import TestCase
import warnings

from dark.genbank import GenomeRanges


class TestErrors(TestCase):
    """
    Test incorrect ways of instantitating the GenomeRanges class.
    """

    def testEmptyString(self):
        """
        The GenomeRanges class must raise a ValueError when passed an range
        empty string.
        """
        error = (
            r'^Could not parse GenBank range string ""\. '
            r'Range "" does not end with \]\(\+\) or \]\(-\)\.$'
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "")

    def testUnparseableString(self):
        """
        The GenomeRanges class must raise a ValueError when passed an
        unparseable string.
        """
        error = (
            r'^Could not parse GenBank range string "x"\. '
            r'Range "x" does not end with \]\(\+\) or \]\(-\)\.$'
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "x")

    def testNoComplementIndicator(self):
        """
        The GenomeRanges class must raise a ValueError when passed an a range
        that does not end with ](+] or ](-).
        """
        error = (
            r'^Could not parse GenBank range string "\[33:40\]"\. '
            r'Range "\[33:40\]" does not end with \]\(\+\) or \]\(-\)\.$'
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "[33:40]")

    def testNoLeadingSquareBracket(self):
        """
        The GenomeRanges class must raise a ValueError when passed an a range
        that does not start with [.
        """
        error = (
            r'^Could not parse GenBank range string "33:40\]\(-\)"\. '
            r'Range "33:40\]\(-\)" does not start with "\["\.$'
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "33:40](-)")

    def testJoinWithNoRange(self):
        """
        The GenomeRanges class must raise a ValueError when passed a join
        directive containing no ranges.
        """
        error = (
            r'^Could not parse GenBank range string "join{}"\. '
            r"join{} can only be used with multiple ranges\.$"
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "join{}")

    def testDecreasingOrder(self):
        """
        The GenomeRanges class must raise a ValueError when passed a string in
        which the start, stop values decrease.
        """
        error = (
            r'^Could not parse GenBank range string "\[7:4\]\(\+\)"\. '
            r"Offset values \(7, 4\) cannot decrease\.$"
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "[7:4](+)")

    def testContiguousRanges(self):
        """
        The GenomeRanges class must raise a ValueError if passed two ranges
        that are contiguous and which have the same (-) direction.
        """
        message = "Contiguous GenBank ranges detected: [3:5] followed by [5:7]."
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            GenomeRanges("join{[3:5](-), [5:7](-)}")
            self.assertEqual(1, len(w))
            self.assertEqual(message, str(w[0].message))

    def testContiguousRangesComplementStrand(self):
        """
        The GenomeRanges class must raise a ValueError if passed two ranges
        that are contiguous and which have the same (+) direction.
        """
        message = "Contiguous GenBank ranges detected: [3:5] followed by [5:7]."
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            GenomeRanges("join{[3:5](+), [5:7](+)}")
            self.assertEqual(1, len(w))
            self.assertEqual(message, str(w[0].message))

    def testTwoNakedRanges(self):
        """
        Two ranges with no join qualifier must result in a ValueError.
        """
        error = (
            r"^Could not parse GenBank range string "
            r'"\[3:5\(\+\), 7:9\(-\)\]"\. '
            r"Multiple ranges must be wrapped in join{}\.$"
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "[3:5(+), 7:9(-)]")

    def testOneJoinedRange(self):
        """
        A single range with a join() must result in a ValueError.
        """
        error = (
            r"^Could not parse GenBank range string "
            r'"join{\[3:5\]\(-\)}"\. '
            r"join{} can only be used with multiple ranges\.$"
        )
        self.assertRaisesRegex(ValueError, error, GenomeRanges, "join{[3:5](-)}")


class TestRanges(TestCase):
    """
    Test the ranges attribute of the GenomeRanges class.
    """

    def testOneNakedRangePositive(self):
        """
        A single range on the positive strand must result in the the expected
        ranges value stored.
        """
        gr = GenomeRanges("[3:5](+)")
        self.assertEqual(((3, 5, True),), gr.ranges)

    def testOneNakedRangeNegative(self):
        """
        A single range on the negative strand must result in the the expected
        ranges value stored.
        """
        gr = GenomeRanges("[3:5](-)")
        self.assertEqual(((3, 5, False),), gr.ranges)

    def testTwoJoinedRanges(self):
        """
        Two joined ranges must return the expected result.
        """
        gr = GenomeRanges("join{[3:5](+), [7:9](-)}")
        self.assertEqual(((3, 5, True), (7, 9, False)), gr.ranges)

    def testThreeJoinedRanges(self):
        """
        Three joined ranges must return the expected result.
        """
        gr = GenomeRanges("join{[3:5](+), [7:9](-), [17:19](-)}")
        self.assertEqual(((3, 5, True), (7, 9, False), (17, 19, False)), gr.ranges)

    def testTwoJoinedContiguousRangesMismatchedStrands(self):
        """
        Two joined ranges that are contiguous but not on the same strand must
        return the expected unmerged (two-range) result.
        """
        gr = GenomeRanges("join{[3:5](-), [5:9](+)}")
        self.assertEqual(((3, 5, False), (5, 9, True)), gr.ranges)

    def testTwoJoinedContiguousRangesComplementStrand(self):
        """
        Two joined ranges that are contiguous on the same (complement) strand
        must return the expected (single-range) result.
        """
        message = "Contiguous GenBank ranges detected: [3:5] followed by [5:9]."
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            gr = GenomeRanges("join{[3:5](-), [5:9](-)}")
            self.assertEqual(((3, 9, False),), gr.ranges)
            self.assertEqual(1, len(w))
            self.assertEqual(message, str(w[0].message))

    def testTwoJoinedContiguousRanges(self):
        """
        Two joined ranges that are contiguous on the same strand must return
        the expected (single-range) result.
        """
        message = "Contiguous GenBank ranges detected: [3:5] followed by [5:9]."
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            gr = GenomeRanges("join{[3:5](+), [5:9](+)}")
            self.assertEqual(((3, 9, True),), gr.ranges)
            self.assertEqual(1, len(w))
            self.assertEqual(message, str(w[0].message))

    def testThreeJoinedContiguousRanges(self):
        """
        Three joined ranges that are contiguous must return the expected
        (single-range) result.
        """
        message1 = "Contiguous GenBank ranges detected: [3:5] followed by [5:7]."
        # Note that the second warning message has a range that doesn't
        # correspond to any of the ranges in the GenomeRanges
        # initialization string. That's because by the time the second
        # warning is issued the first two ranges ([3:5] and [5:7]) have
        # been merged into one ([3:7]).
        message2 = "Contiguous GenBank ranges detected: [3:7] followed by [7:9]."
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            gr = GenomeRanges("join{[3:5](-), [5:7](-), [7:9](-)}")
            self.assertEqual(((3, 9, False),), gr.ranges)
            self.assertEqual(2, len(w))
            self.assertEqual(message1, str(w[0].message))
            self.assertEqual(message2, str(w[1].message))

    def testTwoJoinedContiguousRangesInTheMiddleOfFourRanges(self):
        """
        Two joined ranges that are contiguous must return the expected
        result when there are non-contiguous ranges surrounding them.
        """
        message = "Contiguous GenBank ranges detected: [3:5] followed by [5:9]."
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            gr = GenomeRanges("join{[0:2](-), [3:5](-), [5:9](-), [12:15](-)}")
            self.assertEqual(((0, 2, False), (3, 9, False), (12, 15, False)), gr.ranges)
            self.assertEqual(1, len(w))
            self.assertEqual(message, str(w[0].message))


class TestCircular(TestCase):
    """
    Test the GenomeRanges class ranges method.
    """

    def testOneRange(self):
        """
        The circular method must return False when given only a single range
        tuple that is fully contained within the genome.
        """
        self.assertFalse(GenomeRanges("[20:40](+)").circular(100))

    def testOneRangeEndingAtGenomeEnd(self):
        """
        The circular method must return False when given only a single range
        tuple that ends at the end of the genome.
        """
        self.assertFalse(GenomeRanges("[20:40](+)").circular(40))

    def testOneRangeSpanningTheWholeGenome(self):
        """
        The circular method must return False when given only a single range
        tuple that ends at the end of the genome.
        """
        self.assertFalse(GenomeRanges("[0:40](+)").circular(40))

    def testTwoRangesThatAreCircular(self):
        """
        The circular method must return True when given two ranges that span
        the end of the genome.
        """
        self.assertTrue(GenomeRanges("join{[20:40](+), [0:10](+)}").circular(40))

    def testThreeRangesThatAreCircular(self):
        """
        The circular method must return True when given three ranges that span
        the end of the genome.
        """
        self.assertTrue(
            GenomeRanges("join{[20:30](+), [31:40](+), [0:10](+)}").circular(40)
        )

    def testThreeRangesThatAreCircular2(self):
        """
        The circular method must return True when given three ranges that span
        the end of the genome when one of the ranges is the last in the passed
        BioPython GenBank string and the second range (starting at zero) is the
        first.
        """
        self.assertTrue(
            GenomeRanges("join{[0:10](+), [40:60](+), [80:90](+)}").circular(90)
        )


class TestStartInGenome(TestCase):
    """
    Test the GenomeRanges class startInGenome method.
    """

    def testStartOutsideGenomeRange(self):
        """
        If the start location in the DIAMOND match is greater than the length
        of the genome, a ValueError must be raised.
        """
        # The starting offset in the genome is 27 because DIAMOND returns
        # 1-based offsets into the protein subject. So in Python terms the
        # start offset is amino acid 9, or nucleotide 27.
        ranges = GenomeRanges("[10:20](+)")
        error = (
            r"^Starting nucleotide offset 27 not found in protein "
            r"nucleotide ranges \(10, 20\)\.$"
        )
        self.assertRaisesRegex(ValueError, error, ranges.startInGenome, {"sstart": 10})

    def testStartOutsideGenomeRangeTwoRanges(self):
        """
        If the start location in the DIAMOND match is greater than the length
        of the genome when two ranges are given, a ValueError must be raised.
        """
        # The starting offset in the genome is 27 because DIAMOND returns
        # 1-based offsets into the protein subject. So in Python terms the
        # start offset is amino acid 9, or nucleotide 27.
        ranges = GenomeRanges("join{[5:10](+), [12:18](+)}")
        error = (
            r"^Starting nucleotide offset 27 not found in protein "
            r"nucleotide ranges \(5, 10\), \(12, 18\)\.$"
        )
        self.assertRaisesRegex(ValueError, error, ranges.startInGenome, {"sstart": 10})

    def testMatchAtStartOfFirstRange(self):
        """
        If the start location in the DIAMOND match is at the very beginning of
        the first range, the correct offset must be returned.
        """
        # The offset is 5 because the match starts 0 nucleotides into the
        # protein (0 = (1 - 1 ) * 3) and the protein starts at position 5 in
        # the genome. So the match begins at nucleotide 0 + 5 = 5.
        ranges = GenomeRanges("[5:35](+)")
        self.assertEqual(5, ranges.startInGenome({"sstart": 1}))

    def testMatchAtEndOfFirstRange(self):
        """
        If the start location in the DIAMOND match is at the very end of
        the first range, the correct offset must be returned.
        """
        # The offset is 20 because the match starts 15 nucleotides into the
        # protein (15 = (6 - 1 ) * 3) and the protein starts at position 5 in
        # the genome. So the match begins at nucleotide 15 + 5 = 20.
        ranges = GenomeRanges("[5:21](+)")
        self.assertEqual(20, ranges.startInGenome({"sstart": 6}))

    def testInFirstRange(self):
        """
        If the start location in the DIAMOND match is in the first range,
        the correct offset must be returned.
        """
        # The offset is 32 because the match starts 27 nucleotides into the
        # protein (27 = (10 - 1 ) * 3) and the protein starts at position 5 in
        # the genome. So the match begins at nucleotide 27 + 5 = 32.
        ranges = GenomeRanges("[5:35](+)")
        self.assertEqual(32, ranges.startInGenome({"sstart": 10}))

    def testMatchAtStartOfSecondRange(self):
        """
        If the start location in the DIAMOND match is at the very beginning
        of the second range, the correct offset must be returned.
        """
        # The offset is 50 because the match starts 30 nucleotides into the
        # protein (45 = (16 - 1 ) * 3) and the protein has a range of 30
        # nucleotides (35 - 5 = 30) and then a range of 35 nucleotides (85 -
        # 50 = 35). So the match begins 15 (45 - 30 = 15) nucleotides into
        # the second range (which starts at 50), and 15 + 50 = 65.
        ranges = GenomeRanges("join{[5:35](+), [50:85](+)}")
        self.assertEqual(50, ranges.startInGenome({"sstart": 11}))

    def testInSecondRange(self):
        """
        If the start location in the DIAMOND match is in the second range,
        the correct offset must be returned.
        """
        # The offset is 65 because the match starts 45 nucleotides into the
        # protein (45 = (16 - 1 ) * 3) and the protein has a range of 30
        # nucleotides (35 - 5 = 30) and then a range of 35 nucleotides (85 -
        # 50 = 35). So the match begins 15 (45 - 30 = 15) nucleotides into
        # the second range (which starts at 50), and 15 + 50 = 65.
        ranges = GenomeRanges("join{[5:35](+), [50:85](+)}")
        self.assertEqual(65, ranges.startInGenome({"sstart": 16}))

    def testInThirdRange(self):
        """
        If the start location in the DIAMOND match is in the third range,
        the correct offset must be returned.
        """
        # The offset is 2900 because the match starts 1200 nucleotides into
        # the protein (1200 = (401 - 1 ) * 3) and the protein has a range of
        # 100 nucleotides (100 - 0 = 100), then a range of 200 nucleotides
        # (600 - 400 = 200), then a range of 1000 nucleotides. So the match
        # begins 900 (1200 - 300 = 900) nucleotides into the third range
        # (which starts at 2000), and 2000 + 900 = 2900.
        ranges = GenomeRanges("join{[0:100](+), [400:600](+), [2000:3000](+)}")
        self.assertEqual(2900, ranges.startInGenome({"sstart": 401}))


class TestOrientations(TestCase):
    """
    Test the GenomeRanges class orientations method.
    """

    def testAllComplement(self):
        """
        If all ranges are on the complement strand, the orientations method
        must return just True.
        """
        ranges = GenomeRanges("join{[0:100](-), [400:600](-)}")
        self.assertEqual({False}, ranges.orientations())

    def testNoComplement(self):
        """
        If no ranges are on the complement strand, the orientations method
        must return just False.
        """
        ranges = GenomeRanges("join{[0:100](+), [400:600](+)}")
        self.assertEqual({True}, ranges.orientations())

    def testMixed(self):
        """
        If the ranges are on the both the regular and the complement strand,
        the orientations method must return the set {True, False}.
        """
        ranges = GenomeRanges("join{[0:100](+), [400:600](-)}")
        self.assertEqual({True, False}, ranges.orientations())


class TestDistinctRangeCount(TestCase):
    """
    Test the GenomeRanges class distinctRangeCount method.
    """

    def testOneRange(self):
        """
        If there is only one range, 1 must be returned.
        """
        ranges = GenomeRanges("[0:100](+)")
        self.assertEqual(1, ranges.distinctRangeCount(100))

    def testTwoRange(self):
        """
        If there are 2 ranges, 2 must be returned.
        """
        ranges = GenomeRanges("join{[0:100](+), [400:600](+)}")
        self.assertEqual(2, ranges.distinctRangeCount(100))

    def testThreeRange(self):
        """
        If there are 3 ranges, 3 must be returned.
        """
        ranges = GenomeRanges("join{[0:100](+), [400:600](+), [2000:3000](+)}")
        self.assertEqual(3, ranges.distinctRangeCount(100))

    def testThreeRangesCircular(self):
        """
        If there are 3 ranges but the genome is circular, 2 must be returned.
        """
        ranges = GenomeRanges("join{[0:100](+), [400:600](+), [2000:3000](+)}")
        self.assertEqual(2, ranges.distinctRangeCount(3000))

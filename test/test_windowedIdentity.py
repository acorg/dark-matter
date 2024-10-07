from unittest import TestCase

from dark.reads import DNARead
from dark.windowedIdentity import WindowedIdentity


class TestWindowedIdentity(TestCase):
    """
    Test the WindowedIdentity class.
    """

    def testEmptyFASTA(self):
        """
        If no FASTA sequences are passed, a ValueError must be raised.
        """
        wi = WindowedIdentity([])
        error = "^Empty FASTA input!$"
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "xxx", 10, 3, 2)

    def testWindowSizeTooSmall(self):
        """
        If the window size is too small, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = r"^Window size must be positive\. You passed 0\.$"
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "xxx", 0, 3, 2)

    def testMinWindowSizeNegative(self):
        """
        If the minimum window size is less zero, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = r"^Minimum window size cannot be negative\. You passed -1\.$"
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "xxx", 2, -1, 2)

    def testJumpLessThanOne(self):
        """
        If the jump is less than one, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = r"^Jump must be positive\. You passed 0\.$"
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "xxx", 2, 1, 0)

    def testStartOffsetTooSmall(self):
        """
        If the start offset is too small, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = r"^Start offset cannot be negative\. You passed -1\.$"
        self.assertRaisesRegex(
            ValueError, error, wi.getIdentity, "xxx", 2, 3, 2, startOffset=-1
        )

    def testStopOffsetTooSmall(self):
        """
        If the stop offset is too small, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = r"^Stop offset cannot be negative\. You passed -1\.$"
        self.assertRaisesRegex(
            ValueError, error, wi.getIdentity, "xxx", 2, 3, 2, stopOffset=-1
        )

    def testStopOffsetBeforeStartOffset(self):
        """
        If the stop offset is before the start offset, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = (
            r"^Stop offset must be greater than start offset. You passed startOffset = 4 "
            r"and stopOffset = 3\."
        )
        self.assertRaisesRegex(
            ValueError,
            error,
            wi.getIdentity,
            "xxx",
            2,
            3,
            2,
            startOffset=4,
            stopOffset=3,
        )

    def testStopOffsetEqualsStartOffset(self):
        """
        If the stop offset equals the start offset, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = (
            r"^Stop offset must be greater than start offset. You passed startOffset = 3 "
            r"and stopOffset = 3\."
        )
        self.assertRaisesRegex(
            ValueError,
            error,
            wi.getIdentity,
            "xxx",
            2,
            3,
            2,
            startOffset=3,
            stopOffset=3,
        )

    def testEmptyFirstSequence(self):
        """
        If the first FASTA sequence is empty, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", ""),
            ]
        )
        error = r"^The first input sequence has length zero!$"
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "xxx", 10, 3, 2)

    def testNoMatchingReference(self):
        """
        If the reference pattern doesn't match any sequence id, a ValueError must
        be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = (
            r"^No input sequence IDs match the regular expression for the "
            r"reference \('xxx'\)\.$"
        )
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "xxx", 10, 3, 2)

    def testMultipleMatchingReference(self):
        """
        If the reference pattern matches more than one sequence id, a ValueError must
        be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id1", "ACGT"),
                DNARead("id2", "ACGT"),
            ]
        )
        error = (
            r"^2 input sequence ids match the reference pattern! The matching ids are: "
            r"'id1', 'id2'\.$"
        )
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "id", 10, 3, 2)

    def testNothingToCompareAgainst(self):
        """
        If the include pattern doesn't match any sequence ids, a ValueError must
        be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id", "ACGT"),
            ]
        )
        error = r"^No input sequence ids match the include pattern!$"
        self.assertRaisesRegex(
            ValueError, error, wi.getIdentity, "^id$", 10, 3, 2, includeRegex="xxx"
        )

    def testUnequalReadLengths(self):
        """
        If the passed reads are not all the same length, a ValueError must be raised.
        """
        wi = WindowedIdentity(
            [
                DNARead("id1", "A"),
                DNARead("id2", "AA"),
            ]
        )
        error = (
            r"^Sequence number 2, with id 'id2', has length 2 which is not equal to "
            r"the previous input sequence length\(s\) of 1\.$"
        )
        self.assertRaisesRegex(ValueError, error, wi.getIdentity, "^id$", 10, 3, 2)

    def testPerfectMatchOneIncluded(self):
        """
        If the reference sequence matches the single included sequence perfectly,
        we must get back identities that are all 1.0
        """
        ref = DNARead("reference", "ACGT")
        other = DNARead("another", "ACGT")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 10, 3, 2)
        self.assertIs(ref, reference)
        self.assertEqual([0], starts)
        self.assertEqual({other: [1.0]}, identity)

    def testPerfectMatchOneIncludedStartAtOne(self):
        """
        If the reference sequence matches the single included sequence perfectly,
        we must get back identities that are all 1.0 and the starts list must begin
        at the passed value (of one).
        """
        ref = DNARead("reference", "ACGT")
        other = DNARead("another", "ACGT")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 10, 3, 2, startOffset=1)
        self.assertIs(ref, reference)
        self.assertEqual([1], starts)
        self.assertEqual({other: [1.0]}, identity)

    def testPerfectMatchOneIncludedStartAtOneFirstCharDiffers(self):
        """
        If the reference sequence matches the single included sequence perfectly,
        except for its first character, we must get back identities that are all
        1.0 and the starts list must begin at the passed value (of one).
        """
        ref = DNARead("reference", "CCGT")
        other = DNARead("another", "ACGT")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 10, 3, 2, startOffset=1)
        self.assertIs(ref, reference)
        self.assertEqual([1], starts)
        self.assertEqual({other: [1.0]}, identity)

    def testPerfectMatchOneIncludedIgnoreFirstAndLastChars(self):
        """
        If the reference sequence matches the single included sequence perfectly,
        except for its first and last characters, we must get back identities that
        are all 1.0 and the starts list must begin at the passed value (of one).
        """
        ref = DNARead("reference", "CCGC")
        other = DNARead("another", "ACGT")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity(
            "ref", 10, 2, 2, startOffset=1, stopOffset=3
        )
        self.assertIs(ref, reference)
        self.assertEqual([1], starts)
        self.assertEqual({other: [1.0]}, identity)

    def testOneMatchOneMismatchIgnoreFirstAndLastChars(self):
        """
        If the reference sequence matches the single included sequence in the second
        character and not the third and the first and last characters are ignored and
        the window, min window, and jump size are all one, we must get back identities
        of 1.0 and 0.0 and the starts list must have 1 and 2.
        """
        ref = DNARead("reference", "CCCC")
        other = DNARead("another", "ACGT")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity(
            "ref", 1, 1, 1, startOffset=1, stopOffset=3
        )
        self.assertIs(ref, reference)
        self.assertEqual([1, 2], starts)
        self.assertEqual({other: [1.0, 0.0]}, identity)

    def testMismatchOneIncluded(self):
        """
        If the reference sequence mismatches the single included sequence completely,
        we must get back identities that are all 0.0
        """
        ref = DNARead("reference", "TTTT")
        other = DNARead("another", "GGGG")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 10, 3, 2)
        self.assertIs(ref, reference)
        self.assertEqual([0], starts)
        self.assertEqual({other: [0.0]}, identity)

    def testHalfMatchOneIncluded(self):
        """
        If the reference sequence matches the single included sequence in half its
        locations we must get back identities that are all 0.5
        """
        ref = DNARead("reference", "AATT")
        other = DNARead("another", "AAGG")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 10, 3, 2)
        self.assertIs(ref, reference)
        self.assertEqual([0], starts)
        self.assertEqual({other: [0.5]}, identity)

    def testAmbiguousNonStrictMatch(self):
        """
        If the reference sequence only matches the single included sequence via
        an ambiguous code (just one of two nucelotides matching) we must get back
        an identity of 0.5.
        """
        ref = DNARead("reference", "TT")
        other = DNARead("another", "GW")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 10, 1, 2, strict=False)
        self.assertIs(ref, reference)
        self.assertEqual([0], starts)
        self.assertEqual({other: [0.5]}, identity)

    def testOneIncludedWindowSizeTwo(self):
        """
        If the reference sequence matches the single included sequence perfectly
        in its first two positions, 50% in the second two, and not at all in the final
        two (with a window and jump both equal to two), we must get back identities of
        1.0, 0.5, and 0.0
        """
        ref = DNARead("reference", "AATTCC")
        other = DNARead("another", "AATGAA")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 2, 1, 2)
        self.assertIs(ref, reference)
        self.assertEqual([0, 2, 4], starts)
        self.assertEqual({other: [1.0, 0.5, 0.0]}, identity)

    def testOneIncludedWindowSizeTwoJumpOne(self):
        """
        If the reference sequence matches the single included sequence perfectly
        in its first two positions, 50% in the second two, and not at all in the final
        two (with a window size two, a min window of size two, and jump of one), we must
        get back the expected identities.
        """
        ref = DNARead("reference", "AATTCC")
        other = DNARead("another", "AATGAA")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 2, 2, 1)
        self.assertIs(ref, reference)
        self.assertEqual([0, 1, 2, 3, 4], starts)
        self.assertEqual({other: [1.0, 1.0, 0.5, 0.0, 0.0]}, identity)

    def testOneIncludedWindowSizeTwoJumpOneMinwindowOne(self):
        """
        If the reference sequence matches the single included sequence perfectly
        in its first two positions, 50% in the second two, and not at all in the final
        two (with a window size two, a min window of one and a jump of one), we must
        get back the expected identities.
        """
        ref = DNARead("reference", "AATTCC")
        other = DNARead("another", "AATGAA")
        wi = WindowedIdentity([ref, other])

        reference, starts, identity = wi.getIdentity("ref", 2, 1, 1)
        self.assertIs(ref, reference)
        self.assertEqual([0, 1, 2, 3, 4, 5], starts)
        self.assertEqual({other: [1.0, 1.0, 0.5, 0.0, 0.0, 0.0]}, identity)

    def testOneIncludedWindowSizeTwoIgnoreNonMatchingSequence(self):
        """
        If the reference sequence matches the single included sequence perfectly
        in its first two positions, 50% in the second two, and not at all in the final
        two (with a window and jump both equal to two), we must get back identities of
        1.0, 0.5, and 0.0. An additional sequence with a non-matching id must be
        ignored if it does not match the include regex.
        """
        ref = DNARead("reference", "AATTCC")
        other = DNARead("another", "AATGAA")
        ignored = DNARead("dummy", "AATGAA")
        wi = WindowedIdentity([ref, other, ignored])

        reference, starts, identity = wi.getIdentity("ref", 2, 1, 2, includeRegex="her")
        self.assertIs(ref, reference)
        self.assertEqual([0, 2, 4], starts)
        self.assertEqual({other: [1.0, 0.5, 0.0]}, identity)

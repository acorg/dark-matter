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

    def testDarkMatterIssue802(self):
        """
        Trying to reproduce the failure seen in
        https://github.com/acorg/dark-matter/issues/802 with the exact sequence ids
        that caused that issue.
        """
        ids = """
            NC_003977 FJ899792 JN642140 GQ477453 GQ477455 JN642160 JN642163
            JN688710 JN688711 GQ922005 HE974378 KJ470893 KJ470896 KJ470898
            FJ904430 FJ904436 AB033559 AB048701 AB048702 AB188243 AB210818
            AM494716 AY796031 AY902768 DQ315779 X80925 FJ904399 AY721612
            AY741797 AB270543 EU594409 AB109476 AB555496 GQ205377 MG585269
            AKB003 BRE008 BRE026 BRE028 CUN002 I0157 I0161 I0216 I0217 I1344
            KIL044 KRA001 KRA010 MAY017 MLR005 OAI017 ORE002 Petersberg
            RMI002 SED009 SHK001 STR144 ZWE008 AB486012 AY330911 AJ131571
            AY781180 U46935 AJ131567 AF193863 EU155824 AF498266 BNL005 BON020
            CLL005 DER027 GRG036 I0061 I0104 KLE031 Loschbour MIB003 MIB040
            MIB041 MIS002 MKL025 MN2003 MOT001 MPR001 OOH025 PLZ001 SEC006
            SUA002 TRI011 UOO033 UZZ075 UZZ081 WEH008 WEH009 WSN012 I0461
            KVK001 ARM001 ARM0023 ELT005 GGN005 GRG001 HAL15 HOP004 HSL001
            I0100 I0411 I0551 I0795 I0798 I2020 I2031 JAZ001 MUR007 PDA001
            PEN003 UZZ061 VLI016 YUN048 KAP002 MIB011 AB076679 AB116084
            AB453988 AY738142 GQ477499 AY934764 FJ692556 FJ692598 FJ692611
            GQ161813 GQ331046 AB073858 AB033555 AB219429 AP011089 AB073835
            AB287316 AB287318 AB287320 AB287321 DQ463789 DQ463792 AB241117
            DQ993686 AB111946 AB112066 AB112472 DQ089767 X75656 X75665
            AB048704 AF241411 AP011100 AP011102 AP011103 AP011106 AP011108
            FJ899779 AB049609 FJ023669 EU835241 JN315779 SJN001
            Abusir1543 BOO006 BOO008 CHT001 I1321 KBD002 SGR004 AB048705
            CAO009 AY090458 AB116654 FJ657525 AY090455 AY311369 DQ899144
            DQ899146 AB116549 X75663 AF223962 AB166850 AB059660 AB375163
            AY090454 AY090457 JN792922 KUE033 MAP007 PLS004 SJN013 SGR003
            VLI060 AB056513 AB064312 AF405706 HM363593 X75657 X75664 MK5004
            MK5009 PDA003 LEU065 ELT006 HAL05 ISB002 SID005 UZZ099 VTZ011
            I0784 DA45 DA51 NEO597 NEO238 17662 VK19 VK211-D VK223 RISE254
            DA357 DA195 VK408 Sorsum_low_freq VK18 VK161 DA251
            Sorsum_high_freq DA335 KolymaRiver Karsdorf RISE718 DA222 NEO2324
            DA29 VK477 VK211-G NEO721 DA27 DA119 HBV-A1_N4407 HBV-A2_N4879
            HBV-B2_N4214 HBV-C2_N3825 HBV-D1_N4203 HBV-D3_3091 HBV-G
        """.split()

        ref = DNARead("AB562463", "AATTCC")
        others = [DNARead(id_, "AATGAA") for id_ in ids]
        wi = WindowedIdentity([ref] + others)
        regex = "FJ023669|HBV-G|BOO006|SGR004|AB048705|17662|EU835241|RISE718"

        reference, starts, identity = wi.getIdentity(
            "AB562463", 2, 1, 2, includeRegex=regex
        )

        self.assertIs(ref, reference)
        self.assertEqual([0, 2, 4], starts)
        expected = dict.fromkeys(
            [DNARead(id__, "AATGAA") for id__ in set(regex.split("|"))], [1.0, 0.5, 0.0]
        )
        self.assertEqual(expected, identity)

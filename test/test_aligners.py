from unittest import TestCase

from random import choice, choices

from dark.aligners import edlibAlign, removeUnnecessaryGaps, EDLIB_AMBIGUOUS
from dark.dna import AMBIGUOUS
from dark.reads import DNARead, Reads


class TestEdlibAlign(TestCase):
    """
    Test the edlibAlign function.
    """

    def testAmbiguousCodes(self):
        """
        Test that the edlib ambiguous code tuples are correct.
        """
        # Check a few random tuples (see ../dark/dna.py for more).
        self.assertIn(("A", "M"), EDLIB_AMBIGUOUS)
        self.assertIn(("G", "R"), EDLIB_AMBIGUOUS)
        self.assertNotIn(("R", "G"), EDLIB_AMBIGUOUS)
        self.assertIn(("A", "N"), EDLIB_AMBIGUOUS)

    def testNumberOfAmbiguousCodes(self):
        """
        Test the number of edlib ambiguous code tuples. The 28 comes from
        summing the lengths of all the ambiguous sets of length > 1
        (see ../dark/dna.py).
        """
        expected = sum(
            len(ambiguities)
            for ambiguities in AMBIGUOUS.values()
            if len(ambiguities) > 1
        )

        self.assertEqual(expected, 28)
        self.assertEqual(expected, len(EDLIB_AMBIGUOUS))

    def testEmtpyStrings(self):
        """
        Aligning two empty reads results in an Exception.
        """
        in1 = DNARead("id1", "")
        in2 = DNARead("id2", "")

        error = (
            r"^The object alignResult contains an empty CIGAR string\. "
            r"Users must run align\(\) with task='path'\. Please check "
            r"the input alignResult\.$"
        )
        self.assertRaisesRegex(Exception, error, edlibAlign, Reads([in1, in2]))

    def testTooManySequences(self):
        """
        Passing more than two sequences must result in a ValueError if
        onlyTwoSequences is True.
        """
        in1 = DNARead("id1", "")
        in2 = DNARead("id2", "")
        in3 = DNARead("id3", "")

        error = r"^Passed 1 unexpected extra sequences\.$"
        self.assertRaisesRegex(Exception, error, edlibAlign, Reads([in1, in2, in3]))

    def testFirstReadAlreadyHasAGap(self):
        """
        If the first read has a gap character we should get a ValueError.
        """
        in1 = DNARead("id1", "-A")
        in2 = DNARead("id2", "AA")

        error = r"^Sequence 'id1' contains one or more gap characters '-'\.$"
        self.assertRaisesRegex(Exception, error, edlibAlign, Reads([in1, in2]))

    def testSecondReadAlreadyHasAGap(self):
        """
        If the second read has a gap character we should get a ValueError.
        """
        in1 = DNARead("id1", "AACCT")
        in2 = DNARead("id2", "AA--T")

        error = r"^Sequence 'id2' contains one or more gap characters '-'\.$"
        self.assertRaisesRegex(Exception, error, edlibAlign, Reads([in1, in2]))

    def testIdenticalStringsOfLengthOne(self):
        """
        Aligning identical reads of length one must produce the expected
        result.
        """
        in1 = DNARead("id1", "A")
        in2 = DNARead("id2", "A")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testThreeIdenticalStringsOfLengthOne(self):
        """
        Aligning identical reads of length one must produce the expected
        result, including if we pass three sequences and onlyTwoSequences
        is False.
        """
        in1 = DNARead("id1", "A")
        in2 = DNARead("id2", "A")
        in3 = DNARead("id3", "A")

        out1, out2 = list(edlibAlign(Reads([in1, in2, in3]), onlyTwoSequences=False))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testNonIdenticalStringsOfLengthOne(self):
        """
        Aligning non-identical reads of length one must produce the expected
        result.
        """
        in1 = DNARead("id1", "A")
        in2 = DNARead("id2", "G")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testIdenticalStringsOfLengthTwo(self):
        """
        Aligning identical reads of length two must produce the expected
        result.
        """
        in1 = DNARead("id1", "AG")
        in2 = DNARead("id2", "AG")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testNonIdenticalStringsOfLengthTwo(self):
        """
        Aligning non-identical reads of length two must produce the expected
        result.
        """
        in1 = DNARead("id1", "AT")
        in2 = DNARead("id2", "GC")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testGappedPrefix(self):
        """
        Aligning a read that is a suffix of another must result in gaps being
        placed at the start of the shorter string.
        """
        in1 = DNARead("id1", "ATGCCGTTCA")
        in2 = DNARead("id2", "GCCGTTCA")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual("--GCCGTTCA", out2.sequence)

    def testLongerGappedPrefix(self):
        """
        Aligning a read that is a suffix of another must result in gaps being
        placed at the expected places at the start of the shorter string.
        """
        in1 = DNARead("id1", "ATGCCGTTCA")
        in2 = DNARead("id2", "CGTTCA")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual("---C-GTTCA", out2.sequence)

    def testGappedSuffix(self):
        """
        Aligning a read that is a prefix of another must result in gaps being
        placed at the end of the shorter string.
        """
        in1 = DNARead("id1", "ATGCCGTTCA")
        in2 = DNARead("id2", "ATGCCGTT")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual("ATGCCGTT--", out2.sequence)

    def testLongerGappedSuffix(self):
        """
        Aligning a read that is a prefix of another must result in gaps being
        placed at the expected places at the end of the shorter string.
        """
        in1 = DNARead("id1", "ATGCCGTTCA")
        in2 = DNARead("id2", "ATGCC")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual("ATGCC-----", out2.sequence)

    def testLongerGappedSuffixWithGap(self):
        """
        Aligning a read that is a prefix of another must result in gaps being
        placed at the expected places at the end of the shorter string.
        """
        in1 = DNARead("id1", "ATGCCGTTCA")
        in2 = DNARead("id2", "ATGCCTT")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual("ATGCCGTTCA", out1.sequence)
        self.assertEqual("ATGCC-TT--", out2.sequence)

    def testGapsAtBothEnds(self):
        """
        Aligning two reads that have a common prefix/suffix must result
        in gaps at the start of one and the end of the other.
        """
        in1 = DNARead("id1", "XXXATGCCGTTCA")
        in2 = DNARead("id2", "ATGCCGTTCAYYY")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual("XXXATGCCGTTCA---", out1.sequence)
        self.assertEqual("---ATGCCGTTCAYYY", out2.sequence)

    def testGapsAtBothEndsAndOneInTheMiddle(self):
        """
        Aligning two reads that have a common prefix/suffix and one deletion
        must result in gaps at the start of one and the end of the other,
        and one in the middle.
        """
        in1 = DNARead("id1", "XXXATGCCTTCA")
        in2 = DNARead("id2", "ATGCCGTTCAYYY")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual("XXXATGCC-TTCA---", out1.sequence)
        self.assertEqual("---ATGCCGTTCAYYY", out2.sequence)

    def testGapsAtBothEndsAndTwoInTheMiddle(self):
        """
        Aligning two reads that have a common prefix/suffix and two deletions
        must result in gaps at the start of one and the end of the other,
        and two in the middle.
        """
        in1 = DNARead("id1", "XXXATCCTTCA")
        in2 = DNARead("id2", "ATGCCGTTCAYYY")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual("XXXAT-CC-TTCA---", out1.sequence)
        self.assertEqual("---ATGCCGTTCAYYY", out2.sequence)

    def testOneMiddleGapRemovalNotMinimized(self):
        """
        Test removal of one gap from the middle of each sequence when
        unnecessary gaps are not removed.
        """
        in1 = DNARead("id1", "TAACGCAGA")
        in2 = DNARead("id2", "TAAGTCAGA")

        out1, out2 = list(edlibAlign(Reads([in1, in2]), minimizeGaps=False))

        self.assertEqual("TAACG-CAGA", out1.sequence)
        self.assertEqual("TAA-GTCAGA", out2.sequence)

    def testOneMiddleGapRemovalMinimized(self):
        """
        Test removal of one gap from the middle of each sequence when
        unnecessary gaps are removed.
        """
        in1 = DNARead("id1", "TAACGCAGA")
        in2 = DNARead("id2", "TAAGTCAGA")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testTwoMiddleGapsRemovalNotMinimized(self):
        """
        Test removal of two gaps from the middle of each sequence when the
        post-processing step to remove unnecessary gaps is not done.
        """
        in1 = DNARead("id1", "TAACGCAGAxTAACGCAGA")
        in2 = DNARead("id2", "TAAGTCAGAxTAAGTCAGA")

        out1, out2 = list(edlibAlign(Reads([in1, in2]), minimizeGaps=False))

        self.assertEqual("TAACG-CAGAxTAACG-CAGA", out1.sequence)
        self.assertEqual("TAA-GTCAGAxTAA-GTCAGA", out2.sequence)

    def testTwoMiddleGapsRemovalMinimized(self):
        """
        Test removal of two gaps from the middle of each sequence when
        unnecessary gaps are removed.
        """
        in1 = DNARead("id1", "TAACGCAGAxTAACGCAGA")
        in2 = DNARead("id2", "TAAGTCAGAxTAAGTCAGA")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def test30000Nucleotides(self):
        """
        Aligning two reads of length 30K nucleotides must work (and run
        quickly).
        """
        in1 = DNARead("id1", "A" * 30_000)
        in2 = DNARead("id2", "A" * 30_000)

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testNoWildcardMatch(self):
        """
        If wildcards are not used...
        """
        in1 = DNARead("id1", "ACTG")
        in2 = DNARead("id2", "ATG")

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual("ACTG", out1.sequence)
        self.assertEqual("A-TG", out2.sequence)

    def test500SubstitutionsIn30000Nucleotides(self):
        """
        Test that when 500 random substitutions (though not next to each
        other) are made to a random 30K sequence, that the alignment is as
        expected (i.e., that the sequences are perfectly aligned despite the
        mismatches).
        """
        n = 30_000
        nSubs = 500
        seq1 = choices("ACGT", k=n)
        seq2 = seq1[:]

        # Pick some substitution sites, but not next to one another
        # (otherwise there is a chance gaps are introduced).
        subOffsets = set()
        while len(subOffsets) < nSubs:
            offset = choice(range(n))
            if not (subOffsets & {offset - 1, offset, offset + 1}):
                subOffsets.add(offset)

        # Make some changes in seq2.
        for offset in subOffsets:
            b4 = seq2[offset]
            others = list(set("ACGT") - {seq2[offset]})
            self.assertNotIn(b4, others)
            seq2[offset] = choice(others)

        seq1 = "".join(seq1)
        seq2 = "".join(seq2)

        in1 = DNARead("id1", seq1)
        in2 = DNARead("id2", seq2)

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(len(in1), len(out1))
        self.assertEqual(in1, out1)
        self.assertEqual(seq2, out2.sequence)

    def test500DeletionsIn30000Nucleotides(self):
        """
        Test that when 500 random deletions (though not next to each other) are
        made to a random 30K sequence, that the alignment has gaps at the
        expected sites.
        """
        n = 30_000
        nDels = 500

        # Make a first sequence that doesn't have repeated nucelotides.
        # This makes sure edlib doesn't have an adjacent choice about where
        # to put a gap in a string of identical nucelotides.
        seq1 = []
        while len(seq1) < n:
            nt = choice("ACGT")
            if not seq1 or seq1[-1] != nt:
                seq1.append(nt)

        # Pick some deletion sites, but not next to one another (otherwise
        # it is too hard to know what gaps to test for below).
        delOffsets = set()
        while len(delOffsets) < nDels:
            offset = choice(range(n))
            if not (delOffsets & {offset - 1, offset, offset + 1}):
                delOffsets.add(offset)

        # seq2 will be a copy of seq1, but sans the deleted sites.
        seq2 = []
        for offset in range(n):
            if offset not in delOffsets:
                seq2.append(seq1[offset])

        seq1 = "".join(seq1)
        seq2 = "".join(seq2)

        in1 = DNARead("id1", seq1)
        in2 = DNARead("id2", seq2)

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        # Check that all the sites that were deleted in seq2 now have a gap.
        for offset in range(n):
            if offset in delOffsets:
                # There is (maybe) a small chance the assertion below will
                # fail, in which case the following code can be uncommented
                # to see what is going on. An info line will be printed,
                # with a '*' at offsets where there should be a gap but
                # isn't and a '|' at gap offsets. This is followed by the
                # two aligned sequences, with '-' characters. I'm leaving
                # this code here for now, until we are very sure it doesn't
                # fail due to some pathological selection of substitution
                # sites and nucelotides. If that does happen, reducing
                # nDels or broadening the exclusion zone around each
                # deletion (to more than +/- 1 nt) should greatly reduce
                # the probability of this test failing.
                #
                # if out2.sequence[offset] != '-':
                #     for a in range(n):
                #         print('*' if a == offset else (
                #             '|' if a in delOffsets else ' '), end='')
                #     print()
                #     print(out1.sequence)
                #     print(out2.sequence)
                self.assertEqual("-", out2.sequence[offset])


class TestRemoveUnnecessaryGaps(TestCase):
    """
    Test the removeUnnecessaryGaps function.
    """

    def testEmptyStrings(self):
        """
        No gaps can be removed from two empty strings.
        """
        self.assertEqual(("", ""), removeUnnecessaryGaps("", ""))

    def testStringsWithNoGaps(self):
        """
        No gaps can be removed from two strings with no gaps.
        """
        self.assertEqual(("AG", "AT"), removeUnnecessaryGaps("AG", "AT"))

    def testOneGapOneIdenticalChar(self):
        """
        Identical strings of length two with one gap have the gap removed.
        """
        self.assertEqual(("A", "A"), removeUnnecessaryGaps("A-", "-A"))
        self.assertEqual(("A", "A"), removeUnnecessaryGaps("-A", "A-"))

    def testOneGapOneNonIdenticalChar(self):
        """
        Non-identical strings of length two with one gap have the gap removed.
        """
        self.assertEqual(("A", "B"), removeUnnecessaryGaps("-A", "B-"))
        self.assertEqual(("A", "B"), removeUnnecessaryGaps("A-", "-B"))

    def testFirstGapInFirstString(self):
        """
        A gap can be removed from two strings when the first string is the
        first to have a gap.
        """
        self.assertEqual(
            ("TAAGTCAGA", "TAACGCAGA"),
            removeUnnecessaryGaps("TAA-GTCAGA", "TAACG-CAGA"),
        )

    def testFirstGapInSecondString(self):
        """
        A gap can be removed from two strings when the second string is the
        first to have a gap.
        """
        self.assertEqual(
            ("TAACGCAGA", "TAAGTCAGA"),
            removeUnnecessaryGaps("TAACG-CAGA", "TAA-GTCAGA"),
        )

    def testConsecutiveGaps(self):
        """
        Two consecutive gaps can be removed.
        """
        self.assertEqual(
            ("TAAGTCAGA", "TAACGCAGA"),
            removeUnnecessaryGaps("TAA--GTCAGA", "TAACG--CAGA"),
        )

from unittest import TestCase

from dark.dna import (
    AMBIGUOUS,
    BASES_TO_AMBIGUOUS,
    compareDNAReads,
    matchToString,
    findKozakConsensus,
    FloatBaseCounts,
    sequenceToRegex,
    leastAmbiguous,
    leastAmbiguousFromCounts,
    Bases,
)
from dark.reads import Read, DNARead, DNAKozakRead

# The following are the letters that used to be on
# from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
# IUPACAmbiguousDNA.letters
# But Bio.Alphabet is now deprecated and will be removed.
AMBIGUOUS_DNA_LETTERS = "GATCRYWSMKHBVDN"

AMBIGUOUS_PAIRS = (
    ("AC", "M"),
    ("AG", "R"),
    ("AT", "W"),
    ("GC", "S"),
    ("GT", "K"),
    ("CT", "Y"),
)

AMBIGUOUS_TRIPLES = (("ACG", "V"), ("ACT", "H"), ("AGT", "D"), ("CGT", "B"))


class TestAmbiguousLetters(TestCase):
    """
    Test the ambiguous DNA mapping.
    """

    def testExpectedAmbiguousLetters(self):
        """
        The ambiguous DNA letters must match those given by IUPAC.
        """
        self.assertEqual(sorted(AMBIGUOUS_DNA_LETTERS), sorted(AMBIGUOUS))

    def testExpectedLengthOne(self):
        """
        The unambiguous DNA letters must be in sets by themselves.
        """
        for base in "ACGT":
            self.assertEqual({base}, AMBIGUOUS[base])

    def testExpectedLengthsGreaterThanOne(self):
        """
        The ambiguous DNA letters must be in sets of size greater than one
        and less than 5.
        """
        for base in set(AMBIGUOUS_DNA_LETTERS) - set("ACGT"):
            self.assertTrue(5 > len(AMBIGUOUS[base]) > 1)

    def testAmbiguousLettersAreAllACGT(self):
        """
        The ambiguous DNA letter sets must all be drawn from A, C, G, T.
        """
        for bases in AMBIGUOUS.values():
            self.assertEqual(set(), bases - set("ACGT"))

    def testAmbiguousLettersAreAllSets(self):
        """
        The ambiguous DNA letter sets must be of the set type (this test
        is in lieu of testing that a tuple or list has no repeats).
        """
        for bases in AMBIGUOUS.values():
            self.assertTrue(isinstance(bases, set))


class TestReversedAmbiguousLetters(TestCase):
    """
    Test the reversed ambiguous DNA mapping.
    """

    def testExpectedSingleBase(self):
        """
        The unambiguous bases must map to themselves.
        """
        for base in "ACGT":
            self.assertEqual(base, BASES_TO_AMBIGUOUS[base])

    def testNIsAllBases(self):
        """
        N must map to all bases.
        """
        self.assertEqual("N", BASES_TO_AMBIGUOUS["ACGT"])


class TestCompareDNAReads(TestCase):
    """
    Test the compareDNAReads function.
    """

    def testEmptySequences(self):
        """
        Two empty sequences must compare as expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 0,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", ""), Read("id2", "")),
        )

    def testExactMatch(self):
        """
        Two sequences that match exactly must compare as expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTT"), Read("id2", "ACGTT")),
        )

    def testNoOffsets(self):
        """
        If an empty set of wanted offsets is passed, the result must be empty.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 0,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ATT-T"), Read("id2", "A-GTC"), offsets=set()),
        )

    def testOffsets(self):
        """
        If a set of wanted offsets is passed, the result must be restricted to
        just those offsets.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 1,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(4, "T", "C")],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(
                Read("id1", "ATT-T"), Read("id2", "A-GTC"), offsets=set([0, 4])
            ),
        )

    def testMatchWithAmbiguityButStrict(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the first
        sequence, must compare as expected when we specify matchAmbiguous=False
        to disallow ambiguous matching.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(5, "S", "C")],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(
                Read("id1", "ACGTTS"), Read("id2", "ACGTTC"), matchAmbiguous=False
            ),
        )

    def testMatchWithAmbiguityInFirst(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the first
        sequence, must compare as expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 1,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [(5, "S", "C")],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTTS"), Read("id2", "ACGTTC")),
        )

    def testMatchWithAmbiguityInSecond(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the
        second sequence, must compare as expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 1,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [(5, "C", "S")],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTTC"), Read("id2", "ACGTTS")),
        )

    def testNonMatchingAmbiguityInFirst(self):
        """
        Two sequences that match exactly, apart from one (incompatible)
        ambiguity in the second sequence, must compare as expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(5, "W", "C")],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTTW"), Read("id2", "ACGTTC")),
        )

    def testMatchWithAmbiguityInBoth(self):
        """
        Two sequences that match exactly, apart from one (compatible)
        ambiguity at the same location in the sequence, must compare as
        expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 1,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [(5, "K", "S")],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTTK"), Read("id2", "ACGTTS")),
        )

    def testMatchWithAmbiguityInBothButStrict(self):
        """
        Two sequences that match exactly, apart from one (compatible)
        ambiguity at the same location in the sequence, must compare as
        expected. Strict.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(5, "K", "S")],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(
                Read("id1", "ACGTTK"), Read("id2", "ACGTTS"), matchAmbiguous=False
            ),
        )

    def testMatchWithIncompatibleAmbiguityInBoth(self):
        """
        Two sequences that match exactly, apart from one (incompatible)
        ambiguity at the same location in the sequence, must compare as
        expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(5, "W", "S")],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTTW"), Read("id2", "ACGTTS")),
        )

    def testMatchWithIncompatibleAmbiguityInBothButStrict(self):
        """
        Two sequences that match exactly, apart from one (incompatible)
        ambiguity at the same location in the sequence, must compare as
        expected. Strict.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(5, "W", "S")],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(
                Read("id1", "ACGTTW"), Read("id2", "ACGTTS"), matchAmbiguous=False
            ),
        )

    def testMatchWithIdenticalAmbiguity(self):
        """
        Two sequences that match exactly, including one (identical)
        ambiguity at the same location in the sequence, must compare as
        expected.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 1,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [(5, "N", "N")],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTTN"), Read("id2", "ACGTTN")),
        )

    def testMatchWithIdenticalAmbiguityButStrict(self):
        """
        Two sequences that match exactly, including one (identical)
        ambiguity at the same location in the sequence, must compare as
        expected. Strict.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(5, "N", "N")],
                },
                "read1": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [5],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(
                Read("id1", "ACGTTN"), Read("id2", "ACGTTN"), matchAmbiguous=False
            ),
        )

    def testGapInFirst(self):
        """
        A gap in the first sequence must be dealt with correctly.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 4,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 1,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [2],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "AC-TT"), Read("id2", "ACGTT")),
        )

    def testGapInSecond(self):
        """
        A gap in the second sequence must be dealt with correctly.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 3,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 2,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [1, 2],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTT"), Read("id2", "A--TT")),
        )

    def testNonDefaultGapChars(self):
        """
        We must be able to specify the gap characters.
        """
        for gap in "+$":
            self.assertEqual(
                {
                    "match": {
                        "identicalMatchCount": 3,
                        "ambiguousMatchCount": 0,
                        "gapMismatchCount": 2,
                        "gapGapMismatchCount": 0,
                        "nonGapMismatchCount": 0,
                        "noCoverageCount": 0,
                        "noCoverageNoCoverageCount": 0,
                        "ambiguousMatches": [],
                        "nonGapMismatches": [],
                    },
                    "read1": {
                        "ambiguousOffsets": [],
                        "extraCount": 0,
                        "gapOffsets": [2],
                        "noCoverageOffsets": [],
                    },
                    "read2": {
                        "ambiguousOffsets": [],
                        "extraCount": 0,
                        "gapOffsets": [0],
                        "noCoverageOffsets": [],
                    },
                },
                compareDNAReads(
                    Read("id1", "AC%sTT" % gap),
                    Read("id2", "%sCGTT" % gap),
                    gapChars="+$",
                ),
            )

    def testGapGap(self):
        """
        Coinciding gaps in the sequences must be dealt with correctly.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 2,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 2,
                    "gapGapMismatchCount": 1,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [2, 3],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [1, 2],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "AC--T"), Read("id2", "A--TT")),
        )

    def testGapAmbiguous(self):
        """
        Testing that the ambiguousOffset shows ambiguous characters paired
        with gaps as expected
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 2,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 2,
                    "gapGapMismatchCount": 1,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [1],
                    "extraCount": 0,
                    "gapOffsets": [2, 3],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [3],
                    "extraCount": 0,
                    "gapOffsets": [1, 2],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "AN--T"), Read("id2", "A--NT")),
        )

    def testNoCoverageInFirst(self):
        """
        A no coverage character in the first sequence must be dealt with
        correctly.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 4,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 1,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [2],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(
                Read("id1", "AC?TT"), Read("id2", "ACGTT"), noCoverageChars="?"
            ),
        )

    def testNoCoverageInSecond(self):
        """
        A no coverage in the second sequence must be dealt with correctly.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 3,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 2,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [1, 2],
                },
            },
            compareDNAReads(
                Read("id1", "ACGTT"), Read("id2", "A??TT"), noCoverageChars="?"
            ),
        )

    def testNoCoverageNoCoverage(self):
        """
        Coinciding no coverages in the sequences must be dealt with correctly.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 2,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 2,
                    "noCoverageNoCoverageCount": 1,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [2, 3],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [1, 2],
                },
            },
            compareDNAReads(
                Read("id1", "AC??T"), Read("id2", "A??TT"), noCoverageChars="?"
            ),
        )

    def testNoCoverageAndGaps(self):
        """
        When no coverage and gap characters are both found in the sequences,
        the result must be correct.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 2,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 3,
                    "gapGapMismatchCount": 2,
                    "nonGapMismatchCount": 1,
                    "noCoverageCount": 2,
                    "noCoverageNoCoverageCount": 1,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(10, "T", "A")],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [5, 6, 7],
                    "noCoverageOffsets": [2, 3],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [5, 6, 8, 9],
                    "noCoverageOffsets": [1, 2],
                },
            },
            compareDNAReads(
                Read("id1", "AC??T---CGT"),
                Read("id2", "A??TT--T--A"),
                noCoverageChars="?",
            ),
        )

    def testExtraInFirst(self):
        """
        If the first sequence has extra bases, they must be indicated in the
        extraCount.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 2,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTTCC"), Read("id2", "ACGTT")),
        )

    def testExtraInSecond(self):
        """
        If the second sequence has extra bases, they must be indicated in the
        extraCount.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 2,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTT"), Read("id2", "ACGTTCC")),
        )

    def testExtraAmbiguous(self):
        """
        If the first sequence has extra bases which are ambiguous,they must
        be indicated in the extraCount and in the ambiguousOffset.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 5,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 0,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [],
                },
                "read1": {
                    "ambiguousOffsets": [6],
                    "extraCount": 2,
                    "gapOffsets": [5],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTT-N"), Read("id2", "ACGTT")),
        )

    def testMismatch(self):
        """
        If the sequences have mismatched (non-ambiguous) bases, their count
        must be given correctly in the nonGapMismatchCount.
        """
        self.assertEqual(
            {
                "match": {
                    "identicalMatchCount": 3,
                    "ambiguousMatchCount": 0,
                    "gapMismatchCount": 0,
                    "gapGapMismatchCount": 0,
                    "nonGapMismatchCount": 2,
                    "noCoverageCount": 0,
                    "noCoverageNoCoverageCount": 0,
                    "ambiguousMatches": [],
                    "nonGapMismatches": [(3, "T", "C"), (4, "T", "C")],
                },
                "read1": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
                "read2": {
                    "ambiguousOffsets": [],
                    "extraCount": 0,
                    "gapOffsets": [],
                    "noCoverageOffsets": [],
                },
            },
            compareDNAReads(Read("id1", "ACGTT"), Read("id2", "ACGCC")),
        )


class TestMatchToString(TestCase):
    """
    Test the matchToString function.
    """

    def testMatchWithAmbiguityButStrict(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the first
        sequence, must compare as expected when we specify matchAmbiguous=False
        to disallow ambiguous matching.
        """
        read1 = Read("id1", "ACGTTS")
        read2 = Read("id2", "ACGTTC")
        match = compareDNAReads(read1, read2, matchAmbiguous=False)

        self.assertEqual(
            """\
Exact matches: 5/6 (83.33%)
Ambiguous matches: 0
Mismatches: 1/6 (16.67%)
  Not involving gaps (i.e., conflicts or ambiguities): 1/6 (16.67%)
  Involving a gap in one sequence: 0
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 0
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 6
    Gaps: 0
    No coverage: 0
    Ambiguous: 1/6 (16.67%)
  Id: id2
    Length: 6
    Gaps: 0
    No coverage: 0
    Ambiguous: 0""",
            matchToString(match, read1, read2, matchAmbiguous=False),
        )

    def testMatchWithAmbiguityAndNotStrict(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the first
        sequence, must compare as expected when we specify matchAmbiguous=True
        to allow ambiguous matching.
        """
        read1 = Read("id1", "ACGTTS")
        read2 = Read("id2", "ACGTTC")
        match = compareDNAReads(read1, read2, matchAmbiguous=True)

        self.assertEqual(
            """\
Exact matches: 5/6 (83.33%)
Ambiguous matches: 1/6 (16.67%)
Exact or ambiguous matches: 6/6 (100.00%)
Mismatches: 0
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 0
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 0
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 6
    Gaps: 0
    No coverage: 0
    Ambiguous: 1/6 (16.67%)
  Id: id2
    Length: 6
    Gaps: 0
    No coverage: 0
    Ambiguous: 0""",
            matchToString(match, read1, read2, matchAmbiguous=True),
        )

    def testGapLocations(self):
        """
        Gap locations must be returned correctly.
        """
        read1 = Read("id1", "TTTTTAAAAAAGCGCG")
        read2 = Read("id2", "TTTTT------GCGCG")
        match = compareDNAReads(read1, read2)
        self.assertEqual(
            """\
Exact matches: 10/16 (62.50%)
Ambiguous matches: 0
Mismatches: 6/16 (37.50%)
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 6/16 (37.50%)
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 0
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    No coverage: 0
    Ambiguous: 0
  Id: id2
    Length: 16
    Gaps: 6/16 (37.50%)
    Gap locations (1-based): 6, 7, 8, 9, 10, 11
    No coverage: 0
    Ambiguous: 0""",
            matchToString(match, read1, read2),
        )

    def testAmbiguousMatchLocations(self):
        """
        Ambiguous match locations must be returned correctly.
        """
        read1 = Read("id1", "TWTTTAAAAAAGCGCG")
        read2 = Read("id2", "TTTTT------GCGCK")
        match = compareDNAReads(read1, read2)
        self.assertEqual(
            """\
Exact matches: 8/16 (50.00%)
Ambiguous matches: 2/16 (12.50%)
Exact or ambiguous matches: 10/16 (62.50%)
Mismatches: 6/16 (37.50%)
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 6/16 (37.50%)
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 0
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    No coverage: 0
    Ambiguous: 1/16 (6.25%)
  Id: id2
    Length: 16
    Gaps: 6/16 (37.50%)
    Gap locations (1-based): 6, 7, 8, 9, 10, 11
    No coverage: 0
    Ambiguous: 1/16 (6.25%)
Ambiguous matches:
    2 W T
    16 G K""",
            matchToString(match, read1, read2, includeAmbiguousMatches=True),
        )

    def testNonGapMismatchLocations(self):
        """
        Non-gap mismatch locations must be returned correctly.
        """
        read1 = Read("id1", "TATTTAAAAAAGCGCG")
        read2 = Read("id2", "TTTTT------GCGCC")
        match = compareDNAReads(read1, read2)
        self.maxDiff = None
        self.assertEqual(
            """\
Exact matches: 8/16 (50.00%)
Ambiguous matches: 0
Mismatches: 8/16 (50.00%)
  Not involving gaps (i.e., conflicts): 2/16 (12.50%)
  Involving a gap in one sequence: 6/16 (37.50%)
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 0
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    No coverage: 0
    Ambiguous: 0
  Id: id2
    Length: 16
    Gaps: 6/16 (37.50%)
    Gap locations (1-based): 6, 7, 8, 9, 10, 11
    No coverage: 0
    Ambiguous: 0
Non-gap mismatches:
    2 A T
    16 G C""",
            matchToString(match, read1, read2, includeNonGapMismatches=True),
        )

    def testNoCoverageLocations(self):
        """
        No coverage locations must be returned correctly.
        """
        read1 = Read("id1", "TTTTTAAAAAAGCGCG")
        read2 = Read("id2", "TTTTT??????GCGCG")
        match = compareDNAReads(read1, read2, noCoverageChars="?")
        self.assertEqual(
            """\
Exact matches: 10/16 (62.50%)
Ambiguous matches: 0
Exact matches (ignoring no coverage sites): 10/10 (100.00%)
Ambiguous matches (ignoring no coverage sites): 0
Mismatches: 0
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 0
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 6/16 (37.50%)
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    No coverage: 0
    Ambiguous: 0
  Id: id2
    Length: 16
    Gaps: 0
    No coverage: 6/16 (37.50%)
    No coverage locations (1-based): 6, 7, 8, 9, 10, 11
    Ambiguous: 0""",
            matchToString(match, read1, read2),
        )

    def testExcludeGapLocations(self):
        """
        If gap locations are not wanted, they should not appear in the result
        of a call to matchToString.
        """
        read1 = Read("id1", "TTTTTAAAAAAGCGCG")
        read2 = Read("id2", "TTTTT------GCGCG")
        match = compareDNAReads(read1, read2)
        self.assertEqual(
            """\
Exact matches: 10/16 (62.50%)
Ambiguous matches: 0
Mismatches: 6/16 (37.50%)
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 6/16 (37.50%)
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 0
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    No coverage: 0
    Ambiguous: 0
  Id: id2
    Length: 16
    Gaps: 6/16 (37.50%)
    No coverage: 0
    Ambiguous: 0""",
            matchToString(match, read1, read2, includeGapLocations=False),
        )

    def testExcludeNoCoverageLocations(self):
        """
        If no coverage locations are not wanted, they should not appear in the
        result of a call to matchToString.
        """
        read1 = Read("id1", "TTTTTAAAAAAGCGCG")
        read2 = Read("id2", "TTTTT??????GCGCG")
        match = compareDNAReads(read1, read2, noCoverageChars="?")
        self.assertEqual(
            """\
Exact matches: 10/16 (62.50%)
Ambiguous matches: 0
Exact matches (ignoring no coverage sites): 10/10 (100.00%)
Ambiguous matches (ignoring no coverage sites): 0
Mismatches: 0
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 0
  Involving a gap in both sequences: 0
  Involving no coverage in one sequence: 6/16 (37.50%)
  Involving no coverage in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    No coverage: 0
    Ambiguous: 0
  Id: id2
    Length: 16
    Gaps: 0
    No coverage: 6/16 (37.50%)
    Ambiguous: 0""",
            matchToString(match, read1, read2, includeNoCoverageLocations=False),
        )


class TestFindKozakConsensus(TestCase):
    """
    Test the findKozakConsensus function.
    """

    def testNoSequence(self):
        """
        If no sequence is given, no Kozak sequence should be found.
        """
        read = DNARead("id", "")
        self.assertEqual([], list(findKozakConsensus(read)))

    def testShortSequence(self):
        """
        If a 4 nt long sequence is given, no Kozak sequence should be found.
        """
        read = DNARead("id", "ATTG")
        self.assertEqual([], list(findKozakConsensus(read)))

    def testOneKozakConsensus(self):
        """
        In a given sequence with an exact Kozak consensus sequence, the offset
        and quality percentage should be as expected.
        """
        read = DNARead("id", "ATTGCCGCCATGGGGG")
        expectedKozakRead = DNAKozakRead(read, 3, 13, 100.0)
        (result,) = list(findKozakConsensus(read))
        self.assertEqual(expectedKozakRead, result)

    def testNoKozakConsensus(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead("id", "ATTGCCTCCATGGGGG")
        self.assertEqual([], list(findKozakConsensus(read)))

    def testFindTwoKozakConsensi(self):
        """
        In a given sequence with two Kozak consensuses with different offsets
        and qualities, the output should be as expected.
        """
        read = DNARead("id", "ATTGCCGCCATGGGGGGCCATGG")
        expectedRead1 = DNARead("id", "ATTGCCGCCATGGGGGGCCATGG")
        expectedRead2 = DNARead("id", "ATTGCCGCCATGGGGGGCCATGG")
        expectedKozakRead1 = DNAKozakRead(expectedRead1, 3, 13, 100.0)
        expectedKozakRead2 = DNAKozakRead(expectedRead2, 13, 23, 60.0)

        self.assertEqual(
            [expectedKozakRead1, expectedKozakRead2], list(findKozakConsensus(read))
        )

    def testNoKozakConsensusAtEnd(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead("id", "ATTGCCTCCATGGGGATG")
        self.assertEqual([], list(findKozakConsensus(read)))

    def testKozakConsensusAtEnd(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead("id", "AAAAAAATTGCCGCCATGG")
        expectedKozakRead = DNAKozakRead(read, 9, 19, 100.0)
        (result,) = list(findKozakConsensus(read))
        self.assertEqual(expectedKozakRead, result)

    def testKozakConsensusATGStart(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead("id", "AAAATGGAAAAAAATTGCCGCC")
        self.assertEqual([], list(findKozakConsensus(read)))

    def testKozakConsensusAtEndPart(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead("id", "AAAAAAATTGCCGCCATG")
        self.assertEqual([], list(findKozakConsensus(read)))


class TestFloatBaseCounts(TestCase):
    """
    Test the FloatBaseCounts class.
    """

    def testOneUnambiguousHomogeneous(self):
        """
        If one unambiguous code is passed to FloatBaseCounts, it must be
        considered homogeneous correctly, depending on the passed level.
        """
        counts = FloatBaseCounts("A")
        self.assertTrue(counts.homogeneous(1.0))
        self.assertTrue(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.0))

    def testTwoIdenticalUnambiguoussHomogeneous(self):
        """
        If an unambiguous code is passed to FloatBaseCounts twice, they must be
        considered homogeneous correctly, depending on the passed level.
        """
        counts = FloatBaseCounts("AA")
        self.assertTrue(counts.homogeneous(1.0))
        self.assertTrue(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.0))

    def testTwoDifferentUnambiguoussHomogeneous(self):
        """
        If two different unambiguous codes are passed to FloatBaseCounts,
        they must be considered homogeneous correctly, depending on the passed
        level.
        """
        counts = FloatBaseCounts("AG")
        self.assertFalse(counts.homogeneous(1.0))
        self.assertTrue(counts.homogeneous(0.5))
        self.assertTrue(counts.homogeneous(0.2))

    def testOneAmbiguousHomogeneousM(self):
        """
        If one ambiguous code (M = A, C) is passed to FloatBaseCounts, it must
        be considered homogeneous correctly, depending on the passed level.
        """
        counts = FloatBaseCounts("M")
        self.assertFalse(counts.homogeneous(1.0))
        self.assertFalse(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.0))

    def testOneAmbiguousHomogeneousV(self):
        """
        If one ambiguous code (V = A, C, G) is passed to FloatBaseCounts, it
        must be considered homogeneous correctly, depending on the passed
        level.
        """
        counts = FloatBaseCounts("V")
        self.assertFalse(counts.homogeneous(1.0))
        self.assertFalse(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.2))
        self.assertTrue(counts.homogeneous(0.0))

    def testTwoIdenticalUnambiguoussVariable(self):
        """
        If an unambiguous code is passed to FloatBaseCounts twice, they must be
        considered non-variable.
        """
        counts = FloatBaseCounts("AA")
        self.assertFalse(counts.variable())

    def testTwoDifferentUnambiguoussVariable(self):
        """
        If two different unambiguous codes are passed to FloatBaseCounts,
        they must be considered variable.
        """
        counts = FloatBaseCounts("AT")
        self.assertTrue(counts.variable())

    def testOneUnambiguousOneAmbiguousVariable(self):
        """
        If one unambiguous code and one incompatible ambiguous code are passed
        to FloatBaseCounts, they must be considered confirm variable.
        """
        counts = FloatBaseCounts("AS")
        self.assertTrue(counts.variable(confirm=True))

    def testOneUnambiguousOneAmbiguousNonVariable(self):
        """
        If one unambiguous code and one compatible ambiguous code are passed
        to FloatBaseCounts, they must not be considered confirm variable
        if confirm is True.
        """
        counts = FloatBaseCounts("AM")
        self.assertFalse(counts.variable(confirm=True))

    def testOneUnambiguousOneAmbiguousVariableUnconfirm(self):
        """
        If one unambiguous code and one compatible ambiguous code are passed
        to FloatBaseCounts, they must be considered variable if confirm
        is False.
        """
        counts = FloatBaseCounts("AM")
        self.assertTrue(counts.variable(confirm=False))

    def testOneUnambiguousOneAmbiguousVariableConfirm(self):
        """
        If one unambiguous code and one incompatible ambiguous code are passed
        to FloatBaseCounts, they must be considered variable if confirm
        is True.
        """
        counts = FloatBaseCounts("AY")
        self.assertTrue(counts.variable(confirm=True))

    def testTwoAmbiguousVariableConfirmFalse(self):
        """
        If two compatible but different ambiguous codes are passed to
        FloatBaseCounts, they must be considered variable if confirm
        is False.
        """
        counts = FloatBaseCounts("MR")
        self.assertTrue(counts.variable(confirm=False))

    def testTwoAmbiguousVariableConfirmTrue(self):
        """
        If two compatible but different ambiguous codes are passed to
        FloatBaseCounts, they must not be considered variable if
        confirm is True.
        """
        counts = FloatBaseCounts("MR")
        self.assertFalse(counts.variable(confirm=True))

    def testOneGapVariableUnconfirm(self):
        """
        If one gap and one unambiguous code are passed to FloatBaseCounts,
        they must be considered variable if confirm is False.
        """
        counts = FloatBaseCounts("A-")
        self.assertTrue(counts.variable(confirm=False))

    def testOneGapVariableConfirm(self):
        """
        If one gap and one unambiguous code are passed to FloatBaseCounts,
        they must be considered variable if confirm is True.
        """
        counts = FloatBaseCounts("A-")
        self.assertTrue(counts.variable(confirm=True))

    def testTwoAmbiguousStr(self):
        """
        If two compatible but different ambiguous codes are passed to
        FloatBaseCounts, they must be converted into a string correctly.
        """
        counts = FloatBaseCounts("MR")
        self.assertEqual("A:1.00 C:0.50 G:0.50 (0.500)", str(counts))

    def testTwoAmbiguousStrWithIntegerTotals(self):
        """
        If two 2-way ambiguous codes are passed to
        FloatBaseCounts, they must be converted into a string correctly
        (i.e., with integer counts).
        """
        counts = FloatBaseCounts("MMA")
        self.assertEqual("A:2 C:1 (0.667)", str(counts))

    def testLowerCase(self):
        """
        If two 2-way ambiguous codes are passed to FloatBaseCounts as lower
        case, they must be converted into a string correctly.
        """
        counts = FloatBaseCounts("mm")
        self.assertEqual("A:1 C:1 (0.500)", str(counts))

    def testMixedCase(self):
        """
        If two 2-way ambiguous codes are passed to FloatBaseCounts in mixed
        case, they must be converted into a string correctly.
        """
        counts = FloatBaseCounts("mM")
        self.assertEqual("A:1 C:1 (0.500)", str(counts))

    def testMostFrequentUnambiguous(self):
        """
        If one unambiguous code passed to FloatBaseCounts is most frequent, the
        mostFrequent method must give the expected result.
        """
        counts = FloatBaseCounts("AAACCTGG")
        self.assertEqual({"A"}, counts.mostFrequent())

    def testEquallyFrequent(self):
        """
        If two codes are passed to FloatBaseCounts in equal numbers, the
        mostFrequent method must give the expected result.
        """
        counts = FloatBaseCounts("AAACCCTGG")
        self.assertEqual(set("AC"), counts.mostFrequent())

    def testMostFrequentWithTwoAmbiguousStr(self):
        """
        If two overlapping ambiguous codes are passed to FloatBaseCounts, the
        mostFrequent method must give the expected result.
        """
        counts = FloatBaseCounts("MR")
        self.assertEqual({"A"}, counts.mostFrequent())

    def testHighestFrequencyUnambiguous(self):
        """
        If one unambiguous code passed to FloatBaseCounts is most frequent, the
        highestFrequency method must give the expected result.
        """
        counts = FloatBaseCounts("AAACCTGG")
        self.assertEqual(0.375, counts.highestFrequency())

    def testEqualFrequency(self):
        """
        If two codes are passed to FloatBaseCounts in equal numbers, the
        highestFrequency method must give the expected result.
        """
        counts = FloatBaseCounts("AAAACCCTGG")
        self.assertEqual(0.4, counts.highestFrequency())

    def testHighestFrequencyWithTwoAmbiguousStr(self):
        """
        If two overlapping ambiguous codes are passed to FloatBaseCounts, the
        highestFrequency method must give the expected result.
        """
        counts = FloatBaseCounts("MR")
        self.assertEqual(0.5, counts.highestFrequency())

    def testLength4Unambiguous(self):
        """
        If all unambiguous bases are given to FloatBaseCounts, its length
        must be 4.
        """
        counts = FloatBaseCounts("AAAACCCTGG")
        self.assertEqual(4, len(counts))

    def testLength3Ambiguous(self):
        """
        If overlapping ambiguous codes are passed to FloatBaseCounts, the
        length must give the expected result.
        """
        counts = FloatBaseCounts("MRC")
        self.assertEqual(3, len(counts))


class TestSequenceToRegex(TestCase):
    """
    Test the sequenceToRegex function.
    """

    def testEmpty(self):
        """
        The empty string should result in an empty regex.
        """
        self.assertEqual("", sequenceToRegex(""))

    def testUnambiguous(self):
        """
        An unambiguous string should result in an identical regex.
        """
        self.assertEqual("ACGT", sequenceToRegex("ACGT"))

    def testOneAmbiguous(self):
        """
        One ambiguous characters should result in the expected regex.
        """
        self.assertEqual("[AG]", sequenceToRegex("R"))

    def testTwoAmbiguous(self):
        """
        Two ambiguous characters should result in the expected regex.
        """
        self.assertEqual("[AG][ACG]", sequenceToRegex("RV"))

    def testMixed(self):
        """
        Mixed ambiguous and non-ambiguous characters should result in the
        expected regex.
        """
        self.assertEqual("A[AG]C[ACG]T", sequenceToRegex("ARCVT"))

    def testN(self):
        """
        An 'N' should result in an ACGT regex.
        """
        self.assertEqual("[ACGT]", sequenceToRegex("N"))

    def testQuestionMark(self):
        """
        A '?' should result in an ACGT regex.
        """
        self.assertEqual("[ACGT]", sequenceToRegex("?"))

    def testGap(self):
        """
        A '-' should result in an ACGT regex.
        """
        self.assertEqual("[ACGT]", sequenceToRegex("-"))

    def testWildcard(self):
        """
        An explicit wildcard should result in an ACGT regex.
        """
        self.assertEqual("[ACGT][ACGT][ACGT]", sequenceToRegex("*#!", wildcards="#*!"))

    def testUnknown(self):
        """
        An unknown character should result in a KeyError.
        """
        error = "^'5'$"
        self.assertRaisesRegex(KeyError, error, sequenceToRegex, "5")


class TestLeastAmbiguous(TestCase):
    """
    Test the leastAmbiguous function.
    """

    def testEmpty(self):
        """
        The empty string should result in a KeyError.
        """
        self.assertRaisesRegex(KeyError, "^''$", leastAmbiguous, "")

    def testUnknownNucleotides(self):
        """
        Unknown nucleotides should result in a KeyError.
        """
        self.assertRaisesRegex(KeyError, "^'123'$", leastAmbiguous, "123")

    def testDuplicationsIgnored(self):
        """
        If nucleotides are duplicated, there should be no problem.
        """
        self.assertEqual("A", leastAmbiguous("AAA"))

    def testDuplicationsDifferentCaseIgnored(self):
        """
        If nucleotides are duplicated in different cases, there should be no
        problem.
        """
        self.assertEqual("A", leastAmbiguous("AaaA"))

    def testSingleNucleotides(self):
        """
        Single nucleotides should be returned as themselves.
        """
        for base in "ACGT":
            self.assertEqual(base, leastAmbiguous(base))

    def testTwoNucleotides(self):
        """
        Two nucleotides should be handled correctly.
        """
        for bases, ambiguous in AMBIGUOUS_PAIRS:
            self.assertEqual(ambiguous, leastAmbiguous(bases))

    def testThreeNucleotides(self):
        """
        Three nucleotides should be handled correctly.
        """
        for bases, ambiguous in AMBIGUOUS_TRIPLES:
            self.assertEqual(ambiguous, leastAmbiguous(bases))

    def testFourNucleotides(self):
        """
        All four nucleotides should be handled correctly.
        """
        self.assertEqual("N", leastAmbiguous("AGCT"))

    def testFourNucleotidesOtherOrder(self):
        """
        All four nucleotides should be handled correctly when given in a
        different order.
        """
        self.assertEqual("N", leastAmbiguous("CGTA"))

    def testTuple(self):
        """
        All four nucleotides passed as a tuple should be handled correctly.
        """
        self.assertEqual("N", leastAmbiguous(tuple("AGCT")))

    def testList(self):
        """
        All four nucleotides passed as a list should be handled correctly.
        """
        self.assertEqual("N", leastAmbiguous(list("AGCT")))


class TestLeastAmbiguousFromBases(TestCase):
    """
    Test the leastAmbiguousFromCounts function.
    """

    def testNegativeCount(self):
        """
        A negative count must result in a ValueError.
        """
        self.assertRaisesRegex(
            ValueError,
            r"^Count for base 'A' is negative \(-1\)\.$",
            leastAmbiguousFromCounts,
            {"A": -1},
            0.9,
        )

    def testNegativeThreshold(self):
        """
        A negative threshold must result in a ValueError.
        """
        self.assertRaisesRegex(
            ValueError,
            r"^Threshold cannot be negative \(-0\.9\)\.$",
            leastAmbiguousFromCounts,
            {"A": 3},
            -0.9,
        )

    def testNoCounts(self):
        """
        If an empty dictionary of counts is passed, 'N' must result.
        """
        self.assertEqual("N", leastAmbiguousFromCounts({}, 0.9))

    def testAllCountsZero(self):
        """
        If the counts are all zero, 'N' must result.
        """
        self.assertEqual("N", leastAmbiguousFromCounts({"A": 0, "G": 0}, 0.9))

    def testOneBase(self):
        """
        If there is a single base, it must be returned.
        """
        for base in "ACGT":
            self.assertEqual(base, leastAmbiguousFromCounts({base: 3}, 0.9))

    def testTwoEqual(self):
        """
        If there are two bases with equal counts, the expected ambiguous code
        must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_PAIRS:
            counts = dict.fromkeys(bases, 1)
            self.assertEqual(ambiguous, leastAmbiguousFromCounts(counts, 0.9))

    def testThreeEqual(self):
        """
        If there are three bases with equal counts, the expected ambiguous code
        must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_TRIPLES:
            counts = dict.fromkeys(bases, 1)
            self.assertEqual(ambiguous, leastAmbiguousFromCounts(counts, 0.9))

    def testFourNucleotides(self):
        """
        If all four nucleotides have equal counts, 'N' must be returned.
        """
        counts = dict.fromkeys("AGCT", 1)
        self.assertEqual("N", leastAmbiguousFromCounts(counts, 0.9))

    def testOneOfTwoOverThreshold(self):
        """
        If one base of two is over the threshold, it must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_PAIRS:
            counts = {bases[0]: 3, bases[1]: 1}
            self.assertEqual(bases[0], leastAmbiguousFromCounts(counts, 0.5))

    def testOneOfTwoBelowThreshold(self):
        """
        If one base is dominant but not over the threshold, the ambiguous code
        must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_PAIRS:
            counts = {bases[0]: 3, bases[1]: 1}
            self.assertEqual(ambiguous, leastAmbiguousFromCounts(counts, 0.9))

    def testOneOfThreeOverThreshold(self):
        """
        If one base of three is over the threshold, it must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_TRIPLES:
            counts = {bases[0]: 3, bases[1]: 1, bases[2]: 1}
            self.assertEqual(bases[0], leastAmbiguousFromCounts(counts, 0.5))

    def testTwoOfThreeOverThreshold(self):
        """
        If two bases of three are over the threshold, the code for those two
        must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_TRIPLES:
            counts = {bases[0]: 4, bases[1]: 4, bases[2]: 2}
            self.assertEqual(
                leastAmbiguous(bases[:2]), leastAmbiguousFromCounts(counts, 0.8)
            )

    def testGeneiousExamplesNoTie(self):
        """
        Test the no-tied counts example from
        https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
        """
        counts = {"A": 6, "G": 3, "T": 1}
        self.assertEqual("A", leastAmbiguousFromCounts(counts, 0.4))
        self.assertEqual("R", leastAmbiguousFromCounts(counts, 0.7))
        self.assertEqual("D", leastAmbiguousFromCounts(counts, 0.95))

    def testGeneiousExamplesTie(self):
        """
        Test the tied counts example from
        https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
        """
        counts = {"A": 6, "G": 2, "T": 2}
        self.assertEqual("A", leastAmbiguousFromCounts(counts, 0.4))
        self.assertEqual("D", leastAmbiguousFromCounts(counts, 0.7))
        self.assertEqual("D", leastAmbiguousFromCounts(counts, 0.95))


class TestBases(TestCase):
    """
    Test the Bases class.
    """

    def testEmptyIsFalse(self):
        """
        An empty Bases instance must be considered False.
        """
        self.assertFalse(Bases())

    def testEmptyHasLengthZero(self):
        """
        An empty Bases instance must have length zero.
        """
        self.assertEqual(0, len(Bases()))

    def testLengthOneIsTrue(self):
        """
        If we append a (base, quality) pair, the instance must be considered
        True.
        """
        b = Bases().append("G", 30)
        self.assertTrue(b)

    def testLengthOne(self):
        """
        If we append a (base, quality) pair, the length must be one.
        """
        b = Bases().append("G", 30)
        self.assertEqual(1, len(b))

    def testGetitem(self):
        """
        __getitem__ must work as expected.
        """
        b = Bases().append("G", 30)
        self.assertEqual(30, b["G"])

    def testAdd(self):
        """
        Addition must work as expected.
        """
        b1 = Bases().append("G", 30)
        b2 = Bases().append("G", 20)
        b3 = b1 + b2

        self.assertEqual(1, len(b1))
        self.assertEqual(1, len(b2))
        self.assertEqual(2, len(b3))
        self.assertEqual(50, b3.counts["G"])

    def testIAdd(self):
        """
        In-place addition must work as expected.
        """
        b1 = Bases().append("G", 30)
        b2 = Bases().append("G", 20)
        b2 += b1

        self.assertEqual(1, len(b1))
        self.assertEqual(2, len(b2))
        self.assertEqual(50, b2.counts["G"])

    def testConsensusNoReads(self):
        """
        The consensus method must return the no coverage string if there are
        no reads.
        """
        b = Bases()
        self.assertEqual("Z", b.consensus(0.8, 1, "L", "Z"))

    def testConsensusLowReads(self):
        """
        The consensus method must return the low coverage string if there are
        no reads.
        """
        b = Bases().append("G", 20)
        self.assertEqual("L", b.consensus(0.8, 2, "L", "Z"))

    def testConsensusOneBase(self):
        """
        If there is a single base, it must be returned.
        """
        for base in "ACGT":
            b = Bases().append(base, 20)
            self.assertEqual(base, b.consensus(0.8, 1, "L", "Z"))

    def testConsensusTwoEqual(self):
        """
        If there are two bases with equal counts, the expected ambiguous code
        must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_PAIRS:
            b = Bases().append(bases[0], 20).append(bases[1], 20)
            self.assertEqual(ambiguous, b.consensus(0.8, 1, "L", "Z"))

    def testConsensusThreeEqual(self):
        """
        If there are three bases with equal counts, the expected ambiguous code
        must be returned.
        """
        for bases, ambiguous in AMBIGUOUS_TRIPLES:
            b = Bases().append(bases[0], 20).append(bases[1], 20).append(bases[2], 20)
            self.assertEqual(ambiguous, b.consensus(0.8, 1, "L", "Z"))

    def testConsensusFourNucleotides(self):
        """
        If all four nucleotides have equal counts, 'N' must be returned.
        """
        b = Bases()
        for base in "ACGT":
            b.append(base, 20)
        self.assertEqual("N", b.consensus(0.8, 1, "L", "Z"))

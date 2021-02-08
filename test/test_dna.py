import six
from unittest import TestCase

from dark.dna import (
    AMBIGUOUS, BASES_TO_AMBIGUOUS, compareDNAReads, matchToString,
    findKozakConsensus, FloatBaseCounts, sequenceToRegex)
from dark.reads import Read, DNARead, DNAKozakRead

# The following are the letters that used to be on
# from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
# IUPACAmbiguousDNA.letters
# But Bio.Alphabet is now deprecated and will be removed.
AMBIGUOUS_DNA_LETTERS = 'GATCRYWSMKHBVDN'


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
        for base in 'ACGT':
            self.assertEqual({base}, AMBIGUOUS[base])

    def testExpectedLengthsGreaterThanOne(self):
        """
        The ambiguous DNA letters must be in sets of size greater than one
        and less than 5.
        """
        for base in set(AMBIGUOUS_DNA_LETTERS) - set('ACGT'):
            self.assertTrue(5 > len(AMBIGUOUS[base]) > 1)

    def testAmbiguousLettersAreAllACGT(self):
        """
        The ambiguous DNA letter sets must all be drawn from A, C, G, T.
        """
        for bases in AMBIGUOUS.values():
            self.assertEqual(set(), bases - set('ACGT'))

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
        for base in 'ACGT':
            self.assertEqual(base, BASES_TO_AMBIGUOUS[base])

    def testNIsAllBases(self):
        """
        N must map to all bases.
        """
        self.assertEqual('N', BASES_TO_AMBIGUOUS['ACGT'])


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
                'match': {
                    'identicalMatchCount': 0,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', ''),
                            Read('id2', '')))

    def testExactMatch(self):
        """
        Two sequences that match exactly must compare as expected.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTT'),
                            Read('id2', 'ACGTT')))

    def testNoOffsets(self):
        """
        If an empty set of wanted offsets is passed, the result must be empty.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 0,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ATT-T'),
                            Read('id2', 'A-GTC'), offsets=set()))

    def testOffsets(self):
        """
        If a set of wanted offsets is passed, the result must be restricted to
        just those offsets.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 1,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ATT-T'),
                            Read('id2', 'A-GTC'), offsets=set([0, 4])))

    def testMatchWithAmbiguityButStrict(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the first
        sequence, must compare as expected when we specify matchAmbiguous=False
        to disallow ambiguous matching.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTS'),
                            Read('id2', 'ACGTTC'), matchAmbiguous=False))

    def testMatchWithAmbiguityInFirst(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the first
        sequence, must compare as expected.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 1,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTS'),
                            Read('id2', 'ACGTTC')))

    def testMatchWithAmbiguityInSecond(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the
        second sequence, must compare as expected.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 1,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTC'),
                            Read('id2', 'ACGTTS')))

    def testNonMatchingAmbiguityInFirst(self):
        """
        Two sequences that match exactly, apart from one (incompatible)
        ambiguity in the second sequence, must compare as expected.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTW'),
                            Read('id2', 'ACGTTC')))

    def testMatchWithAmbiguityInBoth(self):
        """
        Two sequences that match exactly, apart from one (compatible)
        ambiguity at the same location in the sequence, must compare as
        expected.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 1,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTK'),
                            Read('id2', 'ACGTTS')))

    def testMatchWithAmbiguityInBothButStrict(self):
        """
        Two sequences that match exactly, apart from one (compatible)
        ambiguity at the same location in the sequence, must compare as
        expected. Strict.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTK'),
                            Read('id2', 'ACGTTS'), matchAmbiguous=False))

    def testMatchWithIncompatibleAmbiguityInBoth(self):
        """
        Two sequences that match exactly, apart from one (incompatible)
        ambiguity at the same location in the sequence, must compare as
        expected.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTW'),
                            Read('id2', 'ACGTTS')))

    def testMatchWithIncompatibleAmbiguityInBothButStrict(self):
        """
        Two sequences that match exactly, apart from one (incompatible)
        ambiguity at the same location in the sequence, must compare as
        expected. Strict.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTW'),
                            Read('id2', 'ACGTTS'), matchAmbiguous=False))

    def testMatchWithIdenticalAmbiguity(self):
        """
        Two sequences that match exactly, including one (identical)
        ambiguity at the same location in the sequence, must compare as
        expected.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 1,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTN'),
                            Read('id2', 'ACGTTN')))

    def testMatchWithIdenticalAmbiguityButStrict(self):
        """
        Two sequences that match exactly, including one (identical)
        ambiguity at the same location in the sequence, must compare as
        expected. Strict.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [5],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTN'),
                            Read('id2', 'ACGTTN'), matchAmbiguous=False))

    def testGapInFirst(self):
        """
        A gap in the first sequence must be dealt with correctly
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 4,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 1,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [2],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'AC-TT'),
                            Read('id2', 'ACGTT')))

    def testGapInSecond(self):
        """
        A gap in the second sequence must be dealt with correctly
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 3,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 2,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [1, 2],
                },
            },
            compareDNAReads(Read('id1', 'ACGTT'),
                            Read('id2', 'A--TT')))

    def testNonDefaultGapChars(self):
        """
        We must be able to specify the gap characters.
        """
        for gap in '+$':
            self.assertEqual(
                {
                    'match': {
                        'identicalMatchCount': 3,
                        'ambiguousMatchCount': 0,
                        'gapMismatchCount': 2,
                        'gapGapMismatchCount': 0,
                        'nonGapMismatchCount': 0,
                    },
                    'read1': {
                        'ambiguousOffsets': [],
                        'extraCount': 0,
                        'gapOffsets': [2],
                    },
                    'read2': {
                        'ambiguousOffsets': [],
                        'extraCount': 0,
                        'gapOffsets': [0],
                    },
                },
                compareDNAReads(Read('id1', 'AC%sTT' % gap),
                                Read('id2', '%sCGTT' % gap), gapChars='+$'))

    def testGapGap(self):
        """
        Coinciding gaps in the sequences must be dealt with correctly
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 2,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 2,
                    'gapGapMismatchCount': 1,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [2, 3],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [1, 2],
                },
            },
            compareDNAReads(Read('id1', 'AC--T'),
                            Read('id2', 'A--TT')))

    def testGapAmbiguous(self):
        """
        Testing that the ambiguousOffset shows ambiguous characters paired
        with gaps as expected
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 2,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 2,
                    'gapGapMismatchCount': 1,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [1],
                    'extraCount': 0,
                    'gapOffsets': [2, 3],
                },
                'read2': {
                    'ambiguousOffsets': [3],
                    'extraCount': 0,
                    'gapOffsets': [1, 2],
                },
            },
            compareDNAReads(Read('id1', 'AN--T'),
                            Read('id2', 'A--NT')))

    def testExtraInFirst(self):
        """
        If the first sequence has extra bases, they must be indicated in the
        extraCount.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 2,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTTCC'),
                            Read('id2', 'ACGTT')))

    def testExtraInSecond(self):
        """
        If the second sequence has extra bases, they must be indicated in the
        extraCount.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 2,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTT'),
                            Read('id2', 'ACGTTCC')))

    def testExtraAmbiguous(self):
        """
        If the first sequence has extra bases which are ambiguous,they must
        be indicated in the extraCount and in the ambiguousOffset.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 5,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'ambiguousOffsets': [6],
                    'extraCount': 2,
                    'gapOffsets': [5],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTT-N'),
                            Read('id2', 'ACGTT')))

    def testMismatch(self):
        """
        If the sequences have mismatched (non-ambiguous) bases, their count
        must be given correctly in the nonGapMismatchCount.
        """
        self.assertEqual(
            {
                'match': {
                    'identicalMatchCount': 3,
                    'ambiguousMatchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 2,
                },
                'read1': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'ambiguousOffsets': [],
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareDNAReads(Read('id1', 'ACGTT'),
                            Read('id2', 'ACGCC')))


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
        read1 = Read('id1', 'ACGTTS')
        read2 = Read('id2', 'ACGTTC')
        match = compareDNAReads(read1, read2, matchAmbiguous=False)

        self.assertEqual(
            '''\
Exact matches: 5/6 (83.33%)
Ambiguous matches: 0
Mismatches: 1/6 (16.67%)
  Not involving gaps (i.e., conflicts or ambiguities): 1/6 (16.67%)
  Involving a gap in one sequence: 0
  Involving a gap in both sequences: 0
  Id: id1
    Length: 6
    Gaps: 0
    Ambiguous: 1/6 (16.67%)
  Id: id2
    Length: 6
    Gaps: 0
    Ambiguous: 0''',
            matchToString(match, read1, read2, matchAmbiguous=False)
        )

    def testMatchWithAmbiguityAndNotStrict(self):
        """
        Two sequences that match exactly, apart from one ambiguity in the first
        sequence, must compare as expected when we specify matchAmbiguous=True
        to allow ambiguous matching.
        """
        read1 = Read('id1', 'ACGTTS')
        read2 = Read('id2', 'ACGTTC')
        match = compareDNAReads(read1, read2, matchAmbiguous=True)

        self.assertEqual(
            '''\
Exact matches: 5/6 (83.33%)
Ambiguous matches: 1/6 (16.67%)
Exact or ambiguous matches: 6/6 (100.00%)
Mismatches: 0
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 0
  Involving a gap in both sequences: 0
  Id: id1
    Length: 6
    Gaps: 0
    Ambiguous: 1/6 (16.67%)
  Id: id2
    Length: 6
    Gaps: 0
    Ambiguous: 0''',
            matchToString(match, read1, read2, matchAmbiguous=True)
        )

    def testGapLocations(self):
        """
        Gap locations must be returned correctly.
        """
        read1 = Read('id1', 'TTTTTAAAAAAGCGCG')
        read2 = Read('id2', 'TTTTT------GCGCG')
        match = compareDNAReads(read1, read2)
        self.maxDiff = None
        self.assertEqual(
            '''\
Exact matches: 10/16 (62.50%)
Ambiguous matches: 0
Mismatches: 6/16 (37.50%)
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 6/16 (37.50%)
  Involving a gap in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    Ambiguous: 0
  Id: id2
    Length: 16
    Gaps: 6/16 (37.50%)
    Gap locations (1-based): 6, 7, 8, 9, 10, 11
    Ambiguous: 0''',
            matchToString(match, read1, read2)
        )

    def testExcludeGapLocations(self):
        """
        If gap locations are not wanted, they should not appear in the result
        of a call to matchToString.
        """
        read1 = Read('id1', 'TTTTTAAAAAAGCGCG')
        read2 = Read('id2', 'TTTTT------GCGCG')
        match = compareDNAReads(read1, read2)
        self.maxDiff = None
        self.assertEqual(
            '''\
Exact matches: 10/16 (62.50%)
Ambiguous matches: 0
Mismatches: 6/16 (37.50%)
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 6/16 (37.50%)
  Involving a gap in both sequences: 0
  Id: id1
    Length: 16
    Gaps: 0
    Ambiguous: 0
  Id: id2
    Length: 16
    Gaps: 6/16 (37.50%)
    Ambiguous: 0''',
            matchToString(match, read1, read2, includeGapLocations=False)
        )


class TestFindKozakConsensus(TestCase):
    """
    Test the findKozakConsensus function.
    """
    def testNoSequence(self):
        """
        If no sequence is given, no Kozak sequence should be found.
        """
        read = DNARead('id', '')
        self.assertEqual([], list(findKozakConsensus(read)))

    def testShortSequence(self):
        """
        If a 4 nt long sequence is given, no Kozak sequence should be found.
        """
        read = DNARead('id', 'ATTG')
        self.assertEqual([], list(findKozakConsensus(read)))

    def testOneKozakConsensus(self):
        """
        In a given sequence with an exact Kozak consensus sequence, the offset
        and quality percentage should be as expected.
        """
        read = DNARead('id', 'ATTGCCGCCATGGGGG')
        expectedKozakRead = DNAKozakRead(read, 3, 13, 100.0)
        (result,) = list(findKozakConsensus(read))
        self.assertEqual(expectedKozakRead, result)

    def testNoKozakConsensus(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead('id', 'ATTGCCTCCATGGGGG')
        self.assertEqual([], list(findKozakConsensus(read)))

    def testFindTwoKozakConsensi(self):
        """
        In a given sequence with two Kozak consensuses with different offsets
        and qualities, the output should be as expected.
        """
        read = DNARead('id', 'ATTGCCGCCATGGGGGGCCATGG')
        expectedRead1 = DNARead('id', 'ATTGCCGCCATGGGGGGCCATGG')
        expectedRead2 = DNARead('id', 'ATTGCCGCCATGGGGGGCCATGG')
        expectedKozakRead1 = DNAKozakRead(expectedRead1, 3, 13, 100.0)
        expectedKozakRead2 = DNAKozakRead(expectedRead2, 13, 23, 60.0)

        self.assertEqual([expectedKozakRead1, expectedKozakRead2],
                         list(findKozakConsensus(read)))

    def testNoKozakConsensusAtEnd(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead('id', 'ATTGCCTCCATGGGGATG')
        self.assertEqual([], list(findKozakConsensus(read)))

    def testKozakConsensusAtEnd(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead('id', 'AAAAAAATTGCCGCCATGG')
        expectedKozakRead = DNAKozakRead(read, 9, 19, 100.0)
        (result,) = list(findKozakConsensus(read))
        self.assertEqual(expectedKozakRead, result)

    def testKozakConsensusATGStart(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead('id', 'AAAATGGAAAAAAATTGCCGCC')
        self.assertEqual([], list(findKozakConsensus(read)))

    def testKozakConsensusAtEndPart(self):
        """
        In a given sequence without a Kozak consensus, the output should be
        as expected.
        """
        read = DNARead('id', 'AAAAAAATTGCCGCCATG')
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
        counts = FloatBaseCounts('A')
        self.assertTrue(counts.homogeneous(1.0))
        self.assertTrue(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.0))

    def testTwoIdenticalUnambiguoussHomogeneous(self):
        """
        If an unambiguous code is passed to FloatBaseCounts twice, they must be
        considered homogeneous correctly, depending on the passed level.
        """
        counts = FloatBaseCounts('AA')
        self.assertTrue(counts.homogeneous(1.0))
        self.assertTrue(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.0))

    def testTwoDifferentUnambiguoussHomogeneous(self):
        """
        If two different unambiguous codes are passed to FloatBaseCounts,
        they must be considered homogeneous correctly, depending on the passed
        level.
        """
        counts = FloatBaseCounts('AG')
        self.assertFalse(counts.homogeneous(1.0))
        self.assertTrue(counts.homogeneous(0.5))
        self.assertTrue(counts.homogeneous(0.2))

    def testOneAmbiguousHomogeneousM(self):
        """
        If one ambiguous code (M = A, C) is passed to FloatBaseCounts, it must
        be considered homogeneous correctly, depending on the passed level.
        """
        counts = FloatBaseCounts('M')
        self.assertFalse(counts.homogeneous(1.0))
        self.assertFalse(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.0))

    def testOneAmbiguousHomogeneousV(self):
        """
        If one ambiguous code (V = A, C, G) is passed to FloatBaseCounts, it
        must be considered homogeneous correctly, depending on the passed
        level.
        """
        counts = FloatBaseCounts('V')
        self.assertFalse(counts.homogeneous(1.0))
        self.assertFalse(counts.homogeneous(0.75))
        self.assertTrue(counts.homogeneous(0.2))
        self.assertTrue(counts.homogeneous(0.0))

    def testTwoIdenticalUnambiguoussVariable(self):
        """
        If an unambiguous code is passed to FloatBaseCounts twice, they must be
        considered non-variable.
        """
        counts = FloatBaseCounts('AA')
        self.assertFalse(counts.variable())

    def testTwoDifferentUnambiguoussVariable(self):
        """
        If two different unambiguous codes are passed to FloatBaseCounts,
        they must be considered variable.
        """
        counts = FloatBaseCounts('AT')
        self.assertTrue(counts.variable())

    def testOneUnambiguousOneAmbiguousVariable(self):
        """
        If one unambiguous code and one incompatible ambiguous code are passed
        to FloatBaseCounts, they must be considered confirm variable.
        """
        counts = FloatBaseCounts('AS')
        self.assertTrue(counts.variable(confirm=True))

    def testOneUnambiguousOneAmbiguousNonVariable(self):
        """
        If one unambiguous code and one compatible ambiguous code are passed
        to FloatBaseCounts, they must not be considered confirm variable
        if confirm is True.
        """
        counts = FloatBaseCounts('AM')
        self.assertFalse(counts.variable(confirm=True))

    def testOneUnambiguousOneAmbiguousVariableUnconfirm(self):
        """
        If one unambiguous code and one compatible ambiguous code are passed
        to FloatBaseCounts, they must be considered variable if confirm
        is False.
        """
        counts = FloatBaseCounts('AM')
        self.assertTrue(counts.variable(confirm=False))

    def testOneUnambiguousOneAmbiguousVariableConfirm(self):
        """
        If one unambiguous code and one incompatible ambiguous code are passed
        to FloatBaseCounts, they must be considered variable if confirm
        is True.
        """
        counts = FloatBaseCounts('AY')
        self.assertTrue(counts.variable(confirm=True))

    def testTwoAmbiguousVariableConfirmFalse(self):
        """
        If two compatible but different ambiguous codes are passed to
        FloatBaseCounts, they must be considered variable if confirm
        is False.
        """
        counts = FloatBaseCounts('MR')
        self.assertTrue(counts.variable(confirm=False))

    def testTwoAmbiguousVariableConfirmTrue(self):
        """
        If two compatible but different ambiguous codes are passed to
        FloatBaseCounts, they must not be considered variable if
        confirm is True.
        """
        counts = FloatBaseCounts('MR')
        self.assertFalse(counts.variable(confirm=True))

    def testOneGapVariableUnconfirm(self):
        """
        If one gap and one unambiguous code are passed to FloatBaseCounts,
        they must be considered variable if confirm is False.
        """
        counts = FloatBaseCounts('A-')
        self.assertTrue(counts.variable(confirm=False))

    def testOneGapVariableConfirm(self):
        """
        If one gap and one unambiguous code are passed to FloatBaseCounts,
        they must be considered variable if confirm is True.
        """
        counts = FloatBaseCounts('A-')
        self.assertTrue(counts.variable(confirm=True))

    def testTwoAmbiguousStr(self):
        """
        If two compatible but different ambiguous codes are passed to
        FloatBaseCounts, they must be converted into a string correctly.
        """
        counts = FloatBaseCounts('MR')
        self.assertEqual('A:1.00 C:0.50 G:0.50 (0.500)', str(counts))

    def testTwoAmbiguousStrWithIntegerTotals(self):
        """
        If two 2-way ambiguous codes are passed to
        FloatBaseCounts, they must be converted into a string correctly
        (i.e., with integer counts).
        """
        counts = FloatBaseCounts('MMA')
        self.assertEqual('A:2 C:1 (0.667)', str(counts))

    def testLowerCase(self):
        """
        If two 2-way ambiguous codes are passed to FloatBaseCounts as lower
        case, they must be converted into a string correctly.
        """
        counts = FloatBaseCounts('mm')
        self.assertEqual('A:1 C:1 (0.500)', str(counts))

    def testMixedCase(self):
        """
        If two 2-way ambiguous codes are passed to FloatBaseCounts in mixed
        case, they must be converted into a string correctly.
        """
        counts = FloatBaseCounts('mM')
        self.assertEqual('A:1 C:1 (0.500)', str(counts))

    def testMostFrequentUnambiguous(self):
        """
        If one unambiguous code passed to FloatBaseCounts is most frequent, the
        mostFrequent method must give the expected result.
        """
        counts = FloatBaseCounts('AAACCTGG')
        self.assertEqual({'A'}, counts.mostFrequent())

    def testEquallyFrequent(self):
        """
        If two codes are passed to FloatBaseCounts in equal numbers, the
        mostFrequent method must give the expected result.
        """
        counts = FloatBaseCounts('AAACCCTGG')
        self.assertEqual(set('AC'), counts.mostFrequent())

    def testMostFrequentWithTwoAmbiguousStr(self):
        """
        If two overlapping ambiguous codes are passed to FloatBaseCounts, the
        mostFrequent method must give the expected result.
        """
        counts = FloatBaseCounts('MR')
        self.assertEqual({'A'}, counts.mostFrequent())

    def testHighestFrequencyUnambiguous(self):
        """
        If one unambiguous code passed to FloatBaseCounts is most frequent, the
        highestFrequency method must give the expected result.
        """
        counts = FloatBaseCounts('AAACCTGG')
        self.assertEqual(0.375, counts.highestFrequency())

    def testEqualFrequency(self):
        """
        If two codes are passed to FloatBaseCounts in equal numbers, the
        highestFrequency method must give the expected result.
        """
        counts = FloatBaseCounts('AAAACCCTGG')
        self.assertEqual(0.4, counts.highestFrequency())

    def testHighestFrequencyWithTwoAmbiguousStr(self):
        """
        If two overlapping ambiguous codes are passed to FloatBaseCounts, the
        highestFrequency method must give the expected result.
        """
        counts = FloatBaseCounts('MR')
        self.assertEqual(0.5, counts.highestFrequency())

    def testLength4Unambiguous(self):
        """
        If all unambiguous bases are given to FloatBaseCounts, its length
        must be 4.
        """
        counts = FloatBaseCounts('AAAACCCTGG')
        self.assertEqual(4, len(counts))

    def testLength3Ambiguous(self):
        """
        If overlapping ambiguous codes are passed to FloatBaseCounts, the
        length must give the expected result.
        """
        counts = FloatBaseCounts('MRC')
        self.assertEqual(3, len(counts))


class TestSequenceToRegex(TestCase):
    """
    Test the sequenceToRegex function.
    """
    def testEmpty(self):
        """
        The empty string should result in an empty regex.
        """
        self.assertEqual('', sequenceToRegex(''))

    def testUnambiguous(self):
        """
        An unambiguous string should result in an identical regex.
        """
        self.assertEqual('ACGT', sequenceToRegex('ACGT'))

    def testOneAmbiguous(self):
        """
        One ambiguous characters should result in the expected regex.
        """
        self.assertEqual('[AG]', sequenceToRegex('R'))

    def testTwoAmbiguous(self):
        """
        Two ambiguous characters should result in the expected regex.
        """
        self.assertEqual('[AG][ACG]', sequenceToRegex('RV'))

    def testMixed(self):
        """
        Mixed ambiguous and non-ambiguous characters should result in the
        expected regex.
        """
        self.assertEqual('A[AG]C[ACG]T', sequenceToRegex('ARCVT'))

    def testN(self):
        """
        An 'N' should result in an ACGT regex.
        """
        self.assertEqual('[ACGT]', sequenceToRegex('N'))

    def testQuestionMark(self):
        """
        A '?' should result in an ACGT regex.
        """
        self.assertEqual('[ACGT]', sequenceToRegex('?'))

    def testGap(self):
        """
        A '-' should result in an ACGT regex.
        """
        self.assertEqual('[ACGT]', sequenceToRegex('-'))

    def testWildcard(self):
        """
        An explicit wildcard should result in an ACGT regex.
        """
        self.assertEqual('[ACGT][ACGT][ACGT]',
                         sequenceToRegex('*#!', wildcards='#*!'))

    def testUnknown(self):
        """
        An unknown character should result in a KeyError.
        """
        error = "^'5'$"
        six.assertRaisesRegex(self, KeyError, error, sequenceToRegex, '5')

from unittest import TestCase

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from dark.dna import (
    AMBIGUOUS, BASES_TO_AMBIGUOUS, compareDNAReads, matchToString,
    findKozakConsensus)
from dark.reads import Read, DNARead, DNAKozakRead


class TestAmbiguousLetters(TestCase):
    """
    Test the ambiguous DNA mapping.
    """
    def testExpectedAmbiguousLetters(self):
        """
        The ambiguous DNA letters must match those given by IUPAC.
        """
        self.assertEqual(sorted(IUPACAmbiguousDNA.letters), sorted(AMBIGUOUS))

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
        for base in set(IUPACAmbiguousDNA.letters) - set('ACGT'):
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

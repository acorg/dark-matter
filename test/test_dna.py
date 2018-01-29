from unittest import TestCase

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from dark.dna import AMBIGUOUS, BASES_TO_AMBIGUOUS, compareDNAReads
from dark.reads import Read


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

    def testGapInFirst(self):
        """
        A gap in the first sequence must be dealt with correctly
        """
        for gap in '?-N':
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
                compareDNAReads(Read('id1', 'AC%sTT' % gap),
                                Read('id2', 'ACGTT')))

    def testGapInSecond(self):
        """
        A gap in the second sequence must be dealt with correctly
        """
        for gap in '?-N':
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
                                Read('id2', 'A%s%sTT' % (gap, gap))))

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

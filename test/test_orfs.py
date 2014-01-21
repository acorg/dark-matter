from unittest import TestCase
from Bio.Seq import Seq

from dark.orfs import findCodons


class TestFindCodons(TestCase):
    """Tests of the findCodons helper. """

    def testNoMatches(self):
        """
        When there are no codons in the sequence, returns an empty list.
        """
        seq = Seq('AAAAAA')
        self.assertEqual([], list(findCodons(seq, set(['ATG', 'AGG']))))

    def testMatchAtStart(self):
        """
        Finds a codon at the start of the sequence.
        """
        seq = Seq('ATGAAA')
        self.assertEqual([0], list(findCodons(seq, set(['ATG', 'AGG']))))

    def testMatchAtEnd(self):
        """
        Finds a codon at the end of the sequence.
        """
        seq = Seq('ATGAAA')
        self.assertEqual([3], list(findCodons(seq, set(['AAA']))))

    def testMatchMultiple(self):
        """
        Finds multiple codons in the sequence.
        """
        seq = Seq('ATGAAAGGGCCC')
        self.assertEqual([0, 9], list(findCodons(seq, set(['ATG', 'CCC']))))

    def testDoesNotFindOutOfFrameMatches(self):
        """
        Does not find matching codons that are in a non-zero frame in the
        sequence.
        """
        seq = Seq('TATGAAAGGGCCC')
        self.assertEqual([], list(findCodons(seq, set(['ATG', 'CCC']))))

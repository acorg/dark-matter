from unittest import TestCase

from dark.codonDistance import codonInformation, findDistance
from dark.aa import CODONS


class TestCodonInformation(TestCase):
    """
    Tests for the codonInformation function.
    """
    def testReturnRightCodon(self):
        """
        Must return the right codon.
        """
        result = codonInformation(['GAT', 'GAC'], ['GCC', 'GCA'])
        self.assertEqual([['GAC', 'GCC']], result[1])
        self.assertEqual([['GAT', 'GCC'], ['GAT', 'GCA'], ['GAC', 'GCA']],
                         result[2])


class TestFindDistance(TestCase):
    """
    Tests for the findDistance function.
    """
    def testFindDistanceThree(self):
        """
        findDistance must return the right distance if codons are different.
        """
        distance = findDistance('ACA', 'TGT')
        self.assertEqual(3, distance)

    def testFindDistanceZero(self):
        """
        findDistance must return the right distance if codons are identical.
        """
        distance = findDistance('ACA', 'ACA')
        self.assertEqual(0, distance)


class TestCODONS(TestCase):
    """
    Tests for the CODONS table.
    """
    def testNumberOfKeys(self):
        """
        The table must contain the right number of keys.
        """
        self.assertEqual(20, len(CODONS))

    def testNumberCodons(self):
        """
        The table must contain the right number of codons.
        """
        self.assertEqual(44, sum(len(codons) for codons in CODONS.values()))

    def testCodonLength(self):
        """
        All codons must be three bases long.
        """
        self.assertTrue(len(codon) == 3 for codon in CODONS.values())

    def testCodonContent(self):
        """
        Codons must only contain the letters A, T, G, C.
        """
        for aa, codons in CODONS.items():
            for codon in codons:
                self.assertTrue(all(letter in 'ACGT' for letter in codon))

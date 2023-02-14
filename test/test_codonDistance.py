from unittest import TestCase

from dark.codonDistance import codonInformation, findDistance


class TestCodonInformation(TestCase):
    """
    Tests for the codonInformation function.
    """

    def testReturnRightCodon(self):
        """
        The function must return the right codon distances.
        """
        result = codonInformation(["GAT", "GAC"], ["GCC", "GCA"])
        self.assertEqual([("GAC", "GCC")], result[1])
        self.assertEqual([("GAT", "GCC"), ("GAT", "GCA"), ("GAC", "GCA")], result[2])

    def testReturnRightCodonWithTransitionCount(self):
        """
        The function must return the right codon distances and the correct
        number of transitions.
        """
        result = codonInformation(["GAT", "GAC"], ["GCC", "AGT"], countTransitions=True)
        self.assertEqual([("GAC", "GCC", 0)], result[1])
        self.assertEqual([("GAT", "GCC", 1), ("GAT", "AGT", 2)], result[2])
        self.assertEqual([("GAC", "AGT", 3)], result[3])


class TestFindDistance(TestCase):
    """
    Tests for the findDistance function.
    """

    def testFindDistanceThree(self):
        """
        findDistance must return the right distance if codons are different.
        """
        distance = findDistance("ACA", "TGT")
        self.assertEqual(3, distance)

    def testFindDistanceZero(self):
        """
        findDistance must return the right distance if codons are identical.
        """
        distance = findDistance("ACA", "ACA")
        self.assertEqual(0, distance)

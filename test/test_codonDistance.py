from unittest import TestCase

from dark.codonDistance import codonInformation, findDistance


class TestCodonInformation(TestCase):
    """
    Tests for the codonInformation function.
    """
    def testRaiseExceptionIfNoValidCode(self):
        """
        If amino acids are given in not valid code, exception must be raised.
        """
        pass

    def testReturnRightCodon(self):
        """
        Must return the right codon.
        """
        result = codonInformation('Asp', 'Ala')
        self.assertEqual(result['codons']['Asp'], ['GAT', 'GAC'])
        self.assertEqual(result['codons']['Ala'], ['GCC', 'GCA'])

    def testReturnRightDistance(self):
        """
        codonInformation must return the right distances.
        """
        result = codonInformation('Asp', 'Ala')
        print result
        self.assertEqual(result['distances'][0], [])
        self.assertEqual(result['distances'][1], [['GAC', 'GCC']])
        self.assertEqual(result['distances'][2], [['GAT', 'GCC'],
                                                  ['GAT', 'GCA'],
                                                  ['GAC', 'GCA']])
        self.assertEqual(result['distances'][3], [])


class TestFindDistance(TestCase):
    """
    Tests for the findDistance function.
    """
    def testFindDistanceThree(self):
        """
        findDistance must return the right distance if codons are equal.
        """
        distance = findDistance('ACA', 'TGT')
        self.assertEqual(distance, 3)

    def testFindDistanceZero(self):
        """
        findDistance must return the right distance if codons are different.
        """
        distance = findDistance('ACA', 'ACA')
        self.assertEqual(distance, 0)

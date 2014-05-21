from unittest import TestCase
from dark import mutations


class TestGetAPOBECFrequencies(TestCase):
    """
    Tests the getAPOBECFrequencies function
    """
    def testCorrectResult(self):
        """
        The function should return the right bases.
        s:  AGTCTAAAGTCATGACTGGTCCCCTTG
        q:      ..AA..T....G....G...
        """
        result = []
        subjectTitle = '333031'
        subjectSequence = 'AGTCTAAAGTCATGACTGGTCCCCTTG'
        subject = '%s \t \t \t %s' % (subjectTitle, subjectSequence)
        result.append(subject)

        queryTitle = 'H9L0NJQ01DXGOQ'
        querySequence = '    ..AA..T....G....G...   '
        query = '%s \t %s' % (queryTitle, querySequence)
        result.append(query)

        bases = mutations.getAPOBECFrequencies(result, 'C', 'G', 'cPattern')

        self.assertEqual({'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 1, 'CCA': 0,
                          'CCC': 0, 'CCG': 0, 'CCT': 0, 'GCA': 0, 'GCC': 0,
                          'GCG': 0, 'GCT': 0, 'TCA': 0, 'TCC': 1, 'TCG': 0,
                          'TCT': 0}, bases)

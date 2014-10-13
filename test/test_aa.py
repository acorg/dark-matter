from unittest import TestCase
from dark.aa import (
    PROPERTIES, ACIDIC, ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC,
    HYDROPHOBIC, HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR, TINY)


class TestAAProperties(TestCase):
    """
    Test the AA properties in dark.aa
    """

    def testCorrectNumberOfAAs(self):
        """
        The PROPERTIES dict must have the correct AA keys.
        """
        # From https://en.wikipedia.org/wiki/Amino_acid
        letters = 'ARNDCEQGHILKMFPSTWYV'
        self.assertEqual(sorted(letters), sorted(PROPERTIES.keys()))

    def testProperyValuesDiffer(self):
        """
        All individual property values must be different.
        """
        individualProperties = (
            ACIDIC, ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC,
            HYDROPHOBIC, HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR,
            TINY)
        self.assertEqual(len(individualProperties),
                         len(set(individualProperties)))

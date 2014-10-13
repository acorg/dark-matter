from unittest import TestCase

from dark.aa import (
    PROPERTIES, ACIDIC, ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC,
    HYDROPHOBIC, HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR, TINY)

# From https://en.wikipedia.org/wiki/Amino_acid
_AA_LETTERS = 'ARNDCEQGHILKMFPSTWYV'

_INDIVIDUAL_PROPERTIES = (
    ACIDIC, ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC,
    HYDROPHOBIC, HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR, TINY)


class TestAAProperties(TestCase):
    """
    Test the AA properties in dark.aa
    """

    def testCorrectNumberOfAAs(self):
        """
        The PROPERTIES dict must have the correct AA keys.
        """
        self.assertEqual(sorted(_AA_LETTERS), sorted(PROPERTIES.keys()))

    def testPropertyValuesDiffer(self):
        """
        All individual property values must be different.
        """
        self.assertEqual(len(_INDIVIDUAL_PROPERTIES),
                         len(set(_INDIVIDUAL_PROPERTIES)))

    def testEachAAHasItsCorrectPropertiesAndNoOthers(self):
        """
        Each amino acid must have the properties expected of it, and no others.
        """
        expected = {
            'I': set([ALIPHATIC, HYDROPHOBIC]),
            'L': set([ALIPHATIC, HYDROPHOBIC]),
            'V': set([ALIPHATIC, HYDROPHOBIC, SMALL]),
            'M': set([HYDROPHOBIC, SULPHUR]),
            'F': set([HYDROPHOBIC, AROMATIC]),
            'A': set([HYDROPHOBIC, SMALL, TINY]),
            'C': set([HYDROPHOBIC, SMALL, TINY, SULPHUR]),
            'T': set([HYDROPHOBIC, SMALL, HYDROXYLIC]),
            'Y': set([HYDROPHOBIC, AROMATIC, POLAR]),
            'W': set([HYDROPHOBIC, AROMATIC, POLAR]),
            'H': set([HYDROPHOBIC, AROMATIC, POLAR, BASIC_POSITIVE]),
            'K': set([HYDROPHOBIC, BASIC_POSITIVE, POLAR]),
            'P': set([HYDROPHILIC, SMALL]),
            'G': set([HYDROPHILIC, SMALL, TINY]),
            'S': set([HYDROPHILIC, SMALL, POLAR, HYDROXYLIC]),
            'N': set([HYDROPHILIC, SMALL, POLAR, ACIDIC]),
            'D': set([HYDROPHILIC, SMALL, POLAR, NEGATIVE]),
            'Q': set([HYDROPHILIC, POLAR, ACIDIC]),
            'E': set([HYDROPHILIC, NEGATIVE, ACIDIC]),
            'R': set([HYDROPHILIC, POLAR, BASIC_POSITIVE]),
        }

        # Make sure we have expected properties for everything.
        self.assertEqual(sorted(_AA_LETTERS), sorted(expected.keys()))

        for aa in _AA_LETTERS:
            for prop in _INDIVIDUAL_PROPERTIES:
                if prop in expected[aa]:
                    self.assertTrue(bool(PROPERTIES[aa] & prop))
                else:
                    self.assertFalse(bool(PROPERTIES[aa] & prop))

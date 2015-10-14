from unittest import TestCase

from dark.aa import (
    PROPERTIES, ALL_PROPERTIES, PROPERTY_NAMES, ACIDIC,
    ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC, HYDROPHOBIC,
    HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR, TINY, NAMES,
    NAMES_TO_ABBREV1, ABBREV3, ABBREV3_TO_ABBREV1, CODONS, AA_LETTERS,
    find, AminoAcid)


class TestAALetters(TestCase):
    """
    Test the AA_LETTERS string in dark.aa
    """

    def testExpectedAALetters(self):
        """
        The AA_LETTERS value must be as expected.
        """
        # From https://en.wikipedia.org/wiki/Amino_acid
        expected = sorted('ARNDCEQGHILKMFPSTWYV')
        self.assertEqual(expected, AA_LETTERS)


class TestProperties(TestCase):
    """
    Test the AA properties in dark.aa
    """

    def testCorrectAAs(self):
        """
        The PROPERTIES dict must have the correct AA keys.
        """
        self.assertEqual(AA_LETTERS, sorted(PROPERTIES.keys()))

    def testAllProperties(self):
        """
        The ALL_PROPERTIES tuple must contain all known properties.
        """
        expected = set([
            ACIDIC, ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC,
            HYDROPHOBIC, HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR,
            TINY])
        self.assertEqual(set(ALL_PROPERTIES), expected)

    def testPropertyValuesDiffer(self):
        """
        All individual property values must be different.
        """
        self.assertEqual(len(ALL_PROPERTIES), len(set(ALL_PROPERTIES)))

    def testEachAAHasItsCorrectPropertiesAndNoOthers(self):
        """
        Each amino acid must have the properties expected of it, and no others.
        """
        expected = {
            'A': set([HYDROPHOBIC, SMALL, TINY]),
            'C': set([HYDROPHOBIC, SMALL, TINY, SULPHUR]),
            'D': set([HYDROPHILIC, SMALL, POLAR, NEGATIVE]),
            'E': set([HYDROPHILIC, NEGATIVE, ACIDIC]),
            'F': set([HYDROPHOBIC, AROMATIC]),
            'G': set([HYDROPHILIC, SMALL, TINY]),
            'H': set([HYDROPHOBIC, AROMATIC, POLAR, BASIC_POSITIVE]),
            'I': set([ALIPHATIC, HYDROPHOBIC]),
            'K': set([HYDROPHOBIC, BASIC_POSITIVE, POLAR]),
            'L': set([ALIPHATIC, HYDROPHOBIC]),
            'M': set([HYDROPHOBIC, SULPHUR]),
            'N': set([HYDROPHILIC, SMALL, POLAR, ACIDIC]),
            'P': set([HYDROPHILIC, SMALL]),
            'Q': set([HYDROPHILIC, POLAR, ACIDIC]),
            'R': set([HYDROPHILIC, POLAR, BASIC_POSITIVE]),
            'S': set([HYDROPHILIC, SMALL, POLAR, HYDROXYLIC]),
            'T': set([HYDROPHOBIC, SMALL, HYDROXYLIC]),
            'V': set([ALIPHATIC, HYDROPHOBIC, SMALL]),
            'W': set([HYDROPHOBIC, AROMATIC, POLAR]),
            'Y': set([HYDROPHOBIC, AROMATIC, POLAR]),
        }

        # Make sure our 'expected' dict (above) has properties for everything.
        self.assertEqual(AA_LETTERS, sorted(expected.keys()))

        for aa in AA_LETTERS:
            for prop in ALL_PROPERTIES:
                if prop in expected[aa]:
                    self.assertTrue(bool(PROPERTIES[aa] & prop))
                else:
                    self.assertFalse(bool(PROPERTIES[aa] & prop))

    def testProperyNamesKeys(self):
        """
        The PROPERTY_NAMES dict must have a key for each property.
        """
        self.assertEqual(sorted(ALL_PROPERTIES), sorted(PROPERTY_NAMES.keys()))

    def testProperyNamesValuesDiffer(self):
        """
        The PROPERTY_NAMES dict must have different values for each key.
        """
        self.assertEqual(len(PROPERTY_NAMES),
                         len(set(PROPERTY_NAMES.values())))


class TestNames(TestCase):
    """
    Test the NAMES dict in dark.aa
    """

    def testCorrectAAs(self):
        """
        The NAMES dict must have the correct AA keys.
        """
        self.assertEqual(AA_LETTERS, sorted(NAMES.keys()))

    def testValuesDiffer(self):
        """
        All individual names must be different.
        """
        names = list(NAMES.values())
        self.assertEqual(len(names), len(set(names)))


class TestNamesToAbbrev1(TestCase):
    """
    Test the NAMES_TO_ABBREV1 dict in dark.aa
    """

    def testCorrectAAs(self):
        """
        The NAMES_TO_ABBREV1 dict must have the correct AA values.
        """
        self.assertEqual(AA_LETTERS,
                         sorted(NAMES_TO_ABBREV1.values()))


class TestAbbrev3(TestCase):
    """
    Test the ABBREV3 dict in dark.aa
    """

    def testCorrectAAs(self):
        """
        The ABBREV3 dict must have the correct AA keys.
        """
        self.assertEqual(AA_LETTERS, sorted(ABBREV3.keys()))

    def testValuesDiffer(self):
        """
        All individual names must be different.
        """
        names = list(ABBREV3.values())
        self.assertEqual(len(names), len(set(names)))


class TestAbbrev3ToAbbrev1(TestCase):
    """
    Test the ABBREV3_TO_ABBREV1 dict in dark.aa
    """

    def testCorrectAAs(self):
        """
        The ABBREV3_TO_ABBREV1 dict must have the correct AA values.
        """
        self.assertEqual(AA_LETTERS,
                         sorted(ABBREV3_TO_ABBREV1.values()))


class TestCodons(TestCase):
    """
    Tests for the CODONS dict.
    """

    def testCorrectAAs(self):
        """
        The CODONS dict must have the correct AA keys.
        """
        self.assertEqual(AA_LETTERS, sorted(CODONS.keys()))

    def testNumberCodons(self):
        """
        The table must contain the right number of codons.
        """
        self.assertEqual(44, sum(len(codons) for codons in CODONS.values()))

    def testCodonLength(self):
        """
        All codons must be three bases long.
        """
        for codons in CODONS.values():
            self.assertTrue(all(len(codon) == 3 for codon in codons))

    def testCodonContent(self):
        """
        Codons must only contain the letters A, T, G, C and X for unknown.
        """
        for codons in CODONS.values():
            for codon in codons:
                self.assertTrue(all(letter in 'ACGTX' for letter in codon))

    def testCodonsAreNotAbbrev3s(self):
        """
        No codon can be the same as an amino acid 3-letter abbreviation (or
        else our find function may not be unambiguous in what it returns).
        """
        for codons in CODONS.values():
            self.assertFalse(
                any(codon.title() in ABBREV3_TO_ABBREV1 for codon in codons))


class TestAminoAcid(TestCase):
    """
    Test the AminoAcid class in dark.aa
    """

    def testExpectedAttributes(self):
        """
        An amino acid instance must have the expected attributes.
        """
        properties = {}
        propertyDetails = {}
        codons = []
        aa = AminoAcid('Alanine', 'Ala', 'A', codons, properties,
                       propertyDetails)
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertIs(codons, aa.codons)
        self.assertIs(properties, aa.properties)
        self.assertIs(propertyDetails, aa.propertyDetails)


class TestFind(TestCase):
    """
    Test the find function in dark.aa
    """

    def testFindUnknown(self):
        """
        find must return C{None} when called with an unrecognized value.
        """
        self.assertEqual(None, find('silly'))

    def testFindByAbbrev1(self):
        """
        It must be possible to find an amino acid by its 1-letter abbreviation.
        """
        aa = find('A')
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(['GCC', 'GCA'], aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxyethilation': -0.265160523187,
                'iep': -0.191489361702,
                'polar_req': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)

    def testFindByAbbrev3(self):
        """
        It must be possible to find an amino acid by its 3-letter abbreviation.
        """
        aa = find('Ala')
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(['GCC', 'GCA'], aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxyethilation': -0.265160523187,
                'iep': -0.191489361702,
                'polar_req': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)

    def testFindByCodon(self):
        """
        It must be possible to find an amino acid by a 3-letter codon.
        """
        aa = find('GCC')
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(['GCC', 'GCA'], aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxyethilation': -0.265160523187,
                'iep': -0.191489361702,
                'polar_req': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)

    def testFindByName(self):
        """
        It must be possible to find an amino acid by its name.
        """
        aa = find('Alanine')
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(['GCC', 'GCA'], aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxyethilation': -0.265160523187,
                'iep': -0.191489361702,
                'polar_req': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)

    def testFindByNameCaseIgnored(self):
        """
        It must be possible to find an amino acid by its name when the name is
        given in mixed case.
        """
        aa = find('alaNIne')
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(['GCC', 'GCA'], aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxyethilation': -0.265160523187,
                'iep': -0.191489361702,
                'polar_req': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)

    def testFindByNameCaseIgnoredNameWithSpace(self):
        """
        It must be possible to find an amino acid by its name when the name is
        given in mixed case, including if the name has a space in it.
        """
        aa = find('asparTIC aCId')
        self.assertEqual('Aspartic acid', aa.name)
        self.assertEqual('Asp', aa.abbrev3)
        self.assertEqual('D', aa.abbrev1)
        self.assertEqual(['GAT', 'GAC'], aa.codons)
        self.assertEqual(HYDROPHILIC | SMALL | POLAR | NEGATIVE, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': -0.818181818182,
                'aromaticity': -1.0,
                'composition': 0.00363636363636,
                'hydrogenation': -0.90243902439,
                'hydropathy': -0.777777777778,
                'hydroxyethilation': -0.348394768133,
                'iep': -1.0,
                'polar_req': 1.0,
                'polarity': 1.0,
                'volume': -0.389221556886,
            },
            aa.propertyDetails)

import six
from unittest import TestCase
from itertools import product

from dark.aa import (
    PROPERTIES, ALL_PROPERTIES, PROPERTY_NAMES, ACIDIC,
    ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC, HYDROPHOBIC,
    HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR, TINY, NAMES,
    NAMES_TO_ABBREV1, ABBREV3, ABBREV3_TO_ABBREV1, CODONS, AA_LETTERS,
    find, AminoAcid, clustersForSequence, propertiesForSequence,
    PROPERTY_CLUSTERS, PROPERTY_DETAILS_RAW, START_CODON, STOP_CODONS,
    compareAaReads, matchToString)
from dark.reads import AARead


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
        self.assertEqual(AA_LETTERS, sorted(PROPERTIES))

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
        self.assertEqual(AA_LETTERS, sorted(expected))

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
        self.assertEqual(sorted(ALL_PROPERTIES), sorted(PROPERTY_NAMES))

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
        self.assertEqual(AA_LETTERS, sorted(NAMES))

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
        self.assertEqual(AA_LETTERS, sorted(ABBREV3))

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
        self.assertEqual(AA_LETTERS, sorted(CODONS))

    def testNumberCodons(self):
        """
        The table must contain 61 codons. This is 4^3 - 3. I.e., all possible
        combinations of 3 bases minus the three stop codons.
        """
        self.assertEqual(61, sum(len(codons) for codons in CODONS.values()))

    def testCodonLength(self):
        """
        All codons must be three bases long.
        """
        for codons in CODONS.values():
            self.assertTrue(all(len(codon) == 3 for codon in codons))

    def testCodonContent(self):
        """
        Codons must only contain the letters A, T, G, C.
        """
        for codons in CODONS.values():
            for codon in codons:
                self.assertTrue(all(letter in 'ACGT' for letter in codon))

    def testCodonsAreNotAbbrev3s(self):
        """
        No codon can be the same as an amino acid 3-letter abbreviation (or
        else our find function may not be unambiguous in what it returns).
        """
        for codons in CODONS.values():
            self.assertFalse(
                any(codon.title() in ABBREV3_TO_ABBREV1 for codon in codons))

    def testDistinct(self):
        """
        All nucleotide triples must be distinct.
        """
        seen = set()
        for codons in CODONS.values():
            for codon in codons:
                seen.add(codon)
        self.assertEqual(61, len(seen))

    def testStartCodon(self):
        """
        The start codon must have the expected value.
        """
        self.assertEqual('ATG', START_CODON)

    def testStartCodonIsMethionine(self):
        """
        The start codon is the same as Methionine.
        """
        self.assertEqual((START_CODON,), CODONS['M'])

    def testStopCodons(self):
        """
        The stop codons must have the expected value.
        """
        self.assertEqual(('TAA', 'TAG', 'TGA'), STOP_CODONS)

    def testStopCodonsNotInCodonTable(self):
        """
        The stop codons must not be in the main table.
        """
        for stop in STOP_CODONS:
            for codons in CODONS.values():
                self.assertNotIn(stop, codons)

    def testAllCodonsPresent(self):
        """
        All possible codons must be present, as either coding for an AA
        or as a stop codon.
        """
        combinations = set(
            ''.join(x) for x in product('ACGT', 'ACGT', 'ACGT'))
        for codons in CODONS.values():
            for codon in codons:
                combinations.remove(codon)

        # Just the stop codons should be left.
        self.assertEqual(set(STOP_CODONS), combinations)


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
        propertyClusters = {}
        aa = AminoAcid('Alanine', 'Ala', 'A', codons, properties,
                       propertyDetails, propertyClusters)
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertIs(codons, aa.codons)
        self.assertIs(properties, aa.properties)
        self.assertIs(propertyDetails, aa.propertyDetails)
        self.assertIs(propertyClusters, aa.propertyClusters)


class TestFind(TestCase):
    """
    Test the find function in dark.aa
    """

    def testFindUnknown(self):
        """
        find must return an empty generator when called with an unrecognized
        value.
        """
        self.assertEqual([], list(find('silly')))

    def testFindByAbbrev1(self):
        """
        It must be possible to find an amino acid by its 1-letter abbreviation.
        """
        aa = list(find('A'))[0]
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(('GCA', 'GCC', 'GCG', 'GCT'), aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxythiolation': -0.265160523187,
                'iep': -0.191489361702,
                'polar requirement': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)
        self.assertEqual(
            {
                'aliphaticity': 1,
                'aromaticity': 1,
                'composition': 1,
                'hydrogenation': 1,
                'hydropathy': 3,
                'hydroxythiolation': 2,
                'iep': 2,
                'polar requirement': 2,
                'polarity': 2,
                'volume': 2,
            },
            aa.propertyClusters)

    def testFindByAbbrev3(self):
        """
        It must be possible to find an amino acid by its 3-letter abbreviation.
        """
        aa = list(find('Ala'))[0]
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(('GCA', 'GCC', 'GCG', 'GCT'), aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxythiolation': -0.265160523187,
                'iep': -0.191489361702,
                'polar requirement': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)
        self.assertEqual(
            {
                'aliphaticity': 1,
                'aromaticity': 1,
                'composition': 1,
                'hydrogenation': 1,
                'hydropathy': 3,
                'hydroxythiolation': 2,
                'iep': 2,
                'polar requirement': 2,
                'polarity': 2,
                'volume': 2,
            },
            aa.propertyClusters)

    def testFindByCodon(self):
        """
        It must be possible to find an amino acid by a 3-letter codon.
        """
        aa = list(find('GCC'))[0]
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(('GCA', 'GCC', 'GCG', 'GCT'), aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxythiolation': -0.265160523187,
                'iep': -0.191489361702,
                'polar requirement': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)
        self.assertEqual(
            {
                'aliphaticity': 1,
                'aromaticity': 1,
                'composition': 1,
                'hydrogenation': 1,
                'hydropathy': 3,
                'hydroxythiolation': 2,
                'iep': 2,
                'polar requirement': 2,
                'polarity': 2,
                'volume': 2,
            },
            aa.propertyClusters)

    def testFindByName(self):
        """
        It must be possible to find an amino acid by its name.
        """
        aa = list(find('Alanine'))[0]
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(('GCA', 'GCC', 'GCG', 'GCT'), aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxythiolation': -0.265160523187,
                'iep': -0.191489361702,
                'polar requirement': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)
        self.assertEqual(
            {
                'aliphaticity': 1,
                'aromaticity': 1,
                'composition': 1,
                'hydrogenation': 1,
                'hydropathy': 3,
                'hydroxythiolation': 2,
                'iep': 2,
                'polar requirement': 2,
                'polarity': 2,
                'volume': 2,
            },
            aa.propertyClusters)

    def testFindByNameCaseIgnored(self):
        """
        It must be possible to find an amino acid by its name when the name is
        given in mixed case.
        """
        aa = list(find('alaNIne'))[0]
        self.assertEqual('Alanine', aa.name)
        self.assertEqual('Ala', aa.abbrev3)
        self.assertEqual('A', aa.abbrev1)
        self.assertEqual(('GCA', 'GCC', 'GCG', 'GCT'), aa.codons)
        self.assertEqual(HYDROPHOBIC | SMALL | TINY, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': 0.305785123967,
                'aromaticity': -0.550128534704,
                'composition': -1.0,
                'hydrogenation': 0.8973042362,
                'hydropathy': 0.4,
                'hydroxythiolation': -0.265160523187,
                'iep': -0.191489361702,
                'polar requirement': -0.463414634146,
                'polarity': -0.20987654321,
                'volume': -0.664670658683,
            },
            aa.propertyDetails)
        self.assertEqual(
            {
                'aliphaticity': 1,
                'aromaticity': 1,
                'composition': 1,
                'hydrogenation': 1,
                'hydropathy': 3,
                'hydroxythiolation': 2,
                'iep': 2,
                'polar requirement': 2,
                'polarity': 2,
                'volume': 2,
            },
            aa.propertyClusters)

    def testFindByNameCaseIgnoredNameWithSpace(self):
        """
        It must be possible to find an amino acid by its name when the name is
        given in mixed case, including if the name has a space in it.
        """
        aa = list(find('asparTIC aCId'))[0]
        self.assertEqual('Aspartic acid', aa.name)
        self.assertEqual('Asp', aa.abbrev3)
        self.assertEqual('D', aa.abbrev1)
        self.assertEqual(('GAC', 'GAT'), aa.codons)
        self.assertEqual(HYDROPHILIC | SMALL | POLAR | NEGATIVE, aa.properties)
        self.assertEqual(
            {
                'aliphaticity': -0.818181818182,
                'aromaticity': -1.0,
                'composition': 0.00363636363636,
                'hydrogenation': -0.90243902439,
                'hydropathy': -0.777777777778,
                'hydroxythiolation': -0.348394768133,
                'iep': -1.0,
                'polar requirement': 1.0,
                'polarity': 1.0,
                'volume': -0.389221556886,
            },
            aa.propertyDetails)
        self.assertEqual(
            {
                'aliphaticity': 1,
                'aromaticity': 1,
                'composition': 2,
                'hydrogenation': 1,
                'hydropathy': 1,
                'hydroxythiolation': 2,
                'iep': 1,
                'polar requirement': 4,
                'polarity': 4,
                'volume': 3,
            },
            aa.propertyClusters)

    def testFindByPartialName(self):
        """
        It must be possible to find an amino acid by a partial name.
        """
        aa = list(find('nine'))[0]
        self.assertEqual('Alanine', aa.name)

    def testFindByPartialNameMixedCase(self):
        """
        It must be possible to find an amino acid by a mixed-case partial name.
        """
        aa = list(find('NiNe'))[0]
        self.assertEqual('Alanine', aa.name)

    def testFindMultipleMatches(self):
        """
        It must be possible to find more than one matching amino acid.
        """
        aa1, aa2 = list(find('acid'))
        self.assertEqual('Aspartic acid', aa1.name)
        self.assertEqual('Glutamic acid', aa2.name)


class TestPropertyClusters(TestCase):
    """
    Tests for the PROPERTY_CLUSTERS dict.
    """
    def testPropertyClusterKeys(self):
        """
        The PROPERTY_CLUSTERS dict must contain the right keys.
        """
        self.assertEqual(AA_LETTERS, sorted(PROPERTY_CLUSTERS))

    def testNumberOfValues(self):
        """
        Each key in PROPERTY_CLUSTERS must have a dict with 10 elements in it
        as the value.
        """
        for propertyNames in PROPERTY_CLUSTERS.values():
            self.assertEqual([
                'aliphaticity', 'aromaticity', 'composition',
                'hydrogenation', 'hydropathy', 'hydroxythiolation',
                'iep', 'polar requirement', 'polarity', 'volume'],
                sorted(propertyNames))

    def testAliphaticity(self):
        """
        Aliphaticity must always be in cluster 1.
        """
        for propertiesDict in PROPERTY_CLUSTERS.values():
            self.assertEqual(1, propertiesDict['aliphaticity'])

    def testPermittedClustersOnly(self):
        """
        All cluster numbers must be in {1, 2, 3, 4, 5}.
        """
        for propertiesDict in PROPERTY_CLUSTERS.values():
            for clusterNr in propertiesDict.values():
                self.assertIn(clusterNr, {1, 2, 3, 4, 5})


class TestPropertyDetailsRaw(TestCase):
    """
    Tests for the PROPERTY_DETAILS_RAW dict.
    """
    def testPropertyDetailsRawKeys(self):
        """
        The PROPERTY_DETAILS_RAW dict must contain the right keys.
        """
        self.assertEqual(AA_LETTERS, sorted(PROPERTY_DETAILS_RAW))

    def testNumberOfValues(self):
        """
        Each key in PROPERTY_DETAILS_RAW must have a dict with 10 elements in
        it as the value.
        """
        for propertyNames in PROPERTY_DETAILS_RAW.values():
            self.assertEqual([
                'aliphaticity', 'aromaticity', 'composition',
                'hydrogenation', 'hydropathy', 'hydroxythiolation',
                'iep', 'polar requirement', 'polarity', 'volume'],
                sorted(propertyNames))

    def testAliphaticity(self):
        """
        Aliphaticity of Alanin must have the right value.
        """
        self.assertEqual(0.239, PROPERTY_DETAILS_RAW['A']['aliphaticity'])

    def testValuesMustBeFloats(self):
        """
        Each key in PROPERTY_DETAILS_RAW must have a dict whose values are
        all floats.
        """
        for properties in PROPERTY_DETAILS_RAW.values():
            self.assertTrue(all(type(v) is float for v in properties.values()))


class TestPropertiesForSequence(TestCase):
    """
    Tests for the propertiesForSequence function in aa.py
    """
    def testUnknownProperty(self):
        """
        A C{ValueError} must be raised if an unknown property name is passed.
        """
        error = 'Unknown property: xxx'
        read = AARead('id', 'RRR')
        six.assertRaisesRegex(self, ValueError, error,
                              propertiesForSequence, read, ['xxx'])

    def testNoProperties(self):
        """
        If no properties are wanted, an empty dict must be returned.
        """
        read = AARead('id', 'RRR')
        self.assertEqual({}, propertiesForSequence(read, []))

    def testOnePropertyEmptySequence(self):
        """
        If one property is wanted but the sequence is empty, a dict with the
        property must be returned, and have an empty list value.
        """
        read = AARead('id', '')
        self.assertEqual(
            {
                'hydropathy': [],
            },
            propertiesForSequence(read, ['hydropathy']))

    def testPropertyNameIsLowercased(self):
        """
        The property name must be lower-cased in the result.
        """
        read = AARead('id', '')
        self.assertTrue('hydropathy' in
                        propertiesForSequence(read, ['HYDROPATHY']))

    def testOneProperty(self):
        """
        If one property is wanted, a dict with the property must be returned,
        and have the expected property values.
        """
        read = AARead('id', 'AI')
        self.assertEqual(
            {
                'hydropathy': [0.4, 1.0],
            },
            propertiesForSequence(read, ['hydropathy']))

    def testTwoProperties(self):
        """
        If two properties are wanted, a dict with the properties must be
        returned, and have the expected property values.
        """
        read = AARead('id', 'AI')
        self.assertEqual(
            {
                'composition': [-1.0, -1.0],
                'hydropathy': [0.4, 1.0],
            },
            propertiesForSequence(read, ['composition', 'hydropathy']))

    def testDuplicatedPropertyName(self):
        """
        If a property name is mentioned more than once, a dict with the
        wanted properties must be returned, and have the expected property
        values.
        """
        read = AARead('id', 'AI')
        self.assertEqual(
            {
                'composition': [-1.0, -1.0],
                'hydropathy': [0.4, 1.0],
            },
            propertiesForSequence(read, ['composition',
                                         'hydropathy', 'hydropathy']))

    def testMissingAminoAcid(self):
        """
        If an unknown amino acid appears in the sequence, its property value
        must be the default (-1.1).
        """
        read = AARead('id', 'XX')
        self.assertEqual(
            {
                'composition': [-1.1, -1.1],
                'hydropathy': [-1.1, -1.1],
            },
            propertiesForSequence(read, ['composition', 'hydropathy']))

    def testMissingAminoAcidWithNonDefaultMissingValue(self):
        """
        If an unknown amino acid appears in the sequence, its property value
        must be the missing AA value that was passed.
        """
        read = AARead('id', 'XX')
        self.assertEqual(
            {
                'composition': [-1.5, -1.5],
                'hydropathy': [-1.5, -1.5],
            },
            propertiesForSequence(read, ['composition', 'hydropathy'],
                                  missingAAValue=-1.5))


class TestClustersForSequence(TestCase):
    """
    Tests for the clustersForSequence function in aa.py
    """
    def testUnknownProperty(self):
        """
        A C{ValueError} must be raised if an unknown property name is passed.
        """
        error = 'Unknown property: xxx'
        read = AARead('id', 'RRR')
        six.assertRaisesRegex(self, ValueError, error,
                              clustersForSequence, read, ['xxx'])

    def testNoProperties(self):
        """
        If no properties are wanted, an empty dict must be returned.
        """
        read = AARead('id', 'RRR')
        self.assertEqual({}, clustersForSequence(read, []))

    def testOnePropertyEmptySequence(self):
        """
        If one property is wanted but the sequence is empty, a dict with the
        property must be returned, and have an empty list value.
        """
        read = AARead('id', '')
        self.assertEqual(
            {
                'hydropathy': [],
            },
            clustersForSequence(read, ['hydropathy']))

    def testPropertyNameIsLowercased(self):
        """
        The property name must be lower-cased in the result.
        """
        read = AARead('id', '')
        self.assertTrue('hydropathy' in
                        clustersForSequence(read, ['HYDROPATHY']))

    def testOneProperty(self):
        """
        If one property is wanted, a dict with the property must be returned,
        and have the expected property cluster number.
        """
        read = AARead('id', 'AI')
        self.assertEqual(
            {
                'hydropathy': [3, 4],
            },
            clustersForSequence(read, ['hydropathy']))

    def testTwoProperties(self):
        """
        If two properties are wanted, a dict with the properties must be
        returned, and have the expected property cluster numbers.
        """
        read = AARead('id', 'AI')
        self.assertEqual(
            {
                'composition': [1, 1],
                'hydropathy': [3, 4],
            },
            clustersForSequence(read, ['composition', 'hydropathy']))

    def testDuplicatedPropertyName(self):
        """
        If a property name is mentioned more than once, a dict with the
        wanted properties must be returned, and have the expected property
        cluster numbers.
        """
        read = AARead('id', 'AI')
        self.assertEqual(
            {
                'composition': [1, 1],
                'hydropathy': [3, 4],
            },
            clustersForSequence(read, ['composition',
                                       'hydropathy', 'hydropathy']))

    def testMissingAminoAcid(self):
        """
        If an unknown amino acid appears in the sequence, its property cluster
        must be the default (0).
        """
        read = AARead('id', 'XX')
        self.assertEqual(
            {
                'composition': [0, 0],
                'hydropathy': [0, 0],
            },
            clustersForSequence(read, ['composition', 'hydropathy']))

    def testMissingAminoAcidWithNonDefaultMissingValue(self):
        """
        If an unknown amino acid appears in the sequence, its property cluster
        must be the missing AA value that was passed.
        """
        read = AARead('id', 'XX')
        self.assertEqual(
            {
                'composition': [10, 10],
                'hydropathy': [10, 10],
            },
            clustersForSequence(read, ['composition', 'hydropathy'],
                                missingAAValue=10))


class TestCompareAaReads(TestCase):
    """
    Test the compareAaReads function.
    """
    def testEmptySequences(self):
        """
        Two empty sequences must compare as expected.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', ''),
                           AARead('id2', '')))

    def testExactMatch(self):
        """
        Two sequences that match exactly must compare as expected.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 5,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', 'GALHN'),
                           AARead('id2', 'GALHN')))

    def testMismatch(self):
        """
        If the sequences have mismatched bases, their count
        must be given correctly in the nonGapMismatchCount.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 3,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 2,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', 'GALYY'),
                           AARead('id2', 'GALHN')))

    def testNoOffsets(self):
        """
        If an empty set of wanted offsets is passed, the result must be empty.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 0,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', 'GAL-N'),
                           AARead('id2', 'G-LHN'), offsets=set()))

    def testOffsets(self):
        """
        If a set of wanted offsets is passed, the result must be restricted to
        just those offsets.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 1,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 1,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', 'GAL-L'),
                           AARead('id2', 'G-LHN'), offsets=set([0, 4])))

    def testGapInFirst(self):
        """
        A gap in the first sequence must be dealt with correctly.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 4,
                    'gapMismatchCount': 1,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [3],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', 'GAL-N'),
                           AARead('id2', 'GALHN')))

    def testGapInSecond(self):
        """
        A gap in the second sequence must be dealt with correctly.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 3,
                    'gapMismatchCount': 2,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [1, 2],
                },
            },
            compareAaReads(AARead('id1', 'GALHN'),
                           AARead('id2', 'G--HN')))

    def testNonDefaultGapChars(self):
        """
        We must be able to specify the gap characters.
        """
        for gap in '+$':
            self.assertEqual(
                {
                    'match': {
                        'matchCount': 3,
                        'gapMismatchCount': 2,
                        'gapGapMismatchCount': 0,
                        'nonGapMismatchCount': 0,
                    },
                    'read1': {
                        'extraCount': 0,
                        'gapOffsets': [2],
                    },
                    'read2': {
                        'extraCount': 0,
                        'gapOffsets': [0],
                    },
                },
                compareAaReads(AARead('id1', 'GA%sHN' % gap),
                               AARead('id2', '%sALHN' % gap), gapChars='+$'))

    def testGapGap(self):
        """
        Coinciding gaps in the sequences must be dealt with correctly
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 2,
                    'gapMismatchCount': 2,
                    'gapGapMismatchCount': 1,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [2, 3],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [1, 2],
                },
            },
            compareAaReads(AARead('id1', 'GA--N'),
                           AARead('id2', 'G--HN')))

    def testExtraInFirst(self):
        """
        If the first sequence has extra bases, they must be indicated in the
        extraCount.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 5,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 2,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', 'GALHNHN'),
                           AARead('id2', 'GALHN')))

    def testExtraInSecond(self):
        """
        If the second sequence has extra bases, they must be indicated in the
        extraCount.
        """
        self.assertEqual(
            {
                'match': {
                    'matchCount': 5,
                    'gapMismatchCount': 0,
                    'gapGapMismatchCount': 0,
                    'nonGapMismatchCount': 0,
                },
                'read1': {
                    'extraCount': 0,
                    'gapOffsets': [],
                },
                'read2': {
                    'extraCount': 2,
                    'gapOffsets': [],
                },
            },
            compareAaReads(AARead('id1', 'GALHN'),
                           AARead('id2', 'GALHNHN')))


class TestMatchToString(TestCase):
    """
    Test the matchToString function.
    """
    def testMismatchAndMatch(self):
        """
        Two sequences containing matches and mismatches must compare as
        expected.
        """
        read1 = AARead('id1', 'GALHNG')
        read2 = AARead('id2', 'GALHNA')
        match = compareAaReads(read1, read2)

        self.assertEqual(
            '''\
Matches: 5/6 (83.33%)
Mismatches: 1/6 (16.67%)
  Not involving gaps (i.e., conflicts): 1/6 (16.67%)
  Involving a gap in one sequence: 0
  Involving a gap in both sequences: 0
  Id: id1
    Length: 6
    Gaps: 0
  Id: id2
    Length: 6
    Gaps: 0''',
            matchToString(match, read1, read2)
        )

    def testGapAndMatch(self):
        """
        Two sequences containing matches and gaps must compare as
        expected.
        """
        read1 = AARead('id1', 'GALHN-')
        read2 = AARead('id2', 'GALHNA')
        match = compareAaReads(read1, read2)

        self.assertEqual(
            '''\
Matches: 5/6 (83.33%)
Mismatches: 1/6 (16.67%)
  Not involving gaps (i.e., conflicts): 0
  Involving a gap in one sequence: 1/6 (16.67%)
  Involving a gap in both sequences: 0
  Id: id1
    Length: 6
    Gaps: 1/6 (16.67%)
    Gap locations (1-based): 6
  Id: id2
    Length: 6
    Gaps: 0''',
            matchToString(match, read1, read2)
        )

from unittest import TestCase
from os.path import dirname, join

import dark
from dark.genomes import GenomeProteinInfo
from dark.civ.proteins import SqliteIndex

TOP = dirname(dirname(dark.__file__))

DB = SqliteIndex(join(TOP, 'test', 'data', 'hbv', 'hbv-proteins.db'))

BAM1 = join(TOP, 'test', 'data', 'hbv', 'query1.bam')
BAM2 = join(TOP, 'test', 'data', 'hbv', 'query2.bam')
BAM3 = join(TOP, 'test', 'data', 'hbv', 'query3.bam')


class TestGenomeProteinInfo(TestCase):
    """
    Test the GenomeProteinInfo class.
    """

    def testLoadReference(self):
        """
        Test that everything is as expected after loading the genome file.
        """
        gpi = GenomeProteinInfo('KJ586809.1', DB, True)
        self.assertEqual(gpi.genome['proteinCount'], len(gpi.proteins))
        nonProteinOffsets = (
            set(range(gpi.genome['length'])) - set(gpi.offsets))
        self.assertEqual(set(), nonProteinOffsets)
        self.assertEqual('KJ586809.1', gpi.genome['accession'])
        self.assertEqual('Hepatitis B virus strain P18, complete genome',
                         gpi.genome['name'])

    def testLoadBAM1(self):
        """
        Test that everything is as expected after loading the BAM1 file.
        """
        gpi = GenomeProteinInfo('KJ586809.1', DB, True)
        gpi.addSAM(BAM1)

        # Genome covered offsets.
        self.assertEqual(200, len(gpi.coveredOffsetCount))
        self.assertEqual(list(range(200)), list(gpi.coveredOffsetCount))
        self.assertEqual([1] * 200, list(gpi.coveredOffsetCount.values()))

        # SAM files.
        self.assertEqual([BAM1], gpi.samFiles)
        self.assertEqual({'query1'}, gpi.readIdsMatchingGenome)

        # Protein accession numbers.
        expected = set(['AJF208%02d.1' % i for i in range(4, 11)])
        self.assertEqual(expected, set(gpi.proteins))

        # Offset 200 is in 4 proteins but is not matched by the query.
        self.assertEqual(set(), gpi.offsets[200]['readIds'])
        expected = set(['AJF208%02d.1' % i for i in range(4, 8)])
        self.assertEqual(expected, gpi.offsets[200]['proteinAccessions'])

        # Offset 0 is in 3 proteins and is matched by the query.
        self.assertEqual({'query1'}, gpi.offsets[0]['readIds'])
        expected = set(['AJF208%02d.1' % i for i in range(4, 7)])
        self.assertEqual(expected, gpi.offsets[0]['proteinAccessions'])

        # Read ids for all proteins.
        self.assertEqual({'query1'}, gpi.readIdsForAllProteins())

        # AJF20804.1 coverage (its ranges are 2306-3221 and 0-1623)
        info = gpi.proteinCoverageInfo('AJF20804.1')
        self.assertEqual(200, info['coveredOffsets'])
        self.assertEqual(200, info['totalBases'])
        self.assertEqual((3221 - 2306) + (1623 - 0), info['ntLength'])
        self.assertEqual({'query1'}, info['readIds'])

    def testLoadBAM12(self):
        """
        Test that everything is as expected after loading the BAM1 and BAM2
        files.
        """
        gpi = GenomeProteinInfo('KJ586809.1', DB, True)
        gpi.addSAM(BAM1)
        gpi.addSAM(BAM2)

        # Genome covered offsets.
        self.assertEqual(300, len(gpi.coveredOffsetCount))
        self.assertEqual(list(range(200)) + list(range(1400, 1500)),
                         list(gpi.coveredOffsetCount))
        self.assertEqual([1] * 300, list(gpi.coveredOffsetCount.values()))

        # SAM files.
        self.assertEqual([BAM1, BAM2], gpi.samFiles)
        self.assertEqual({'query1', 'query2'}, gpi.readIdsMatchingGenome)

        # Protein accession numbers.
        expected = set(['AJF208%02d.1' % i for i in range(4, 11)])
        self.assertEqual(expected, set(gpi.proteins))

        # Offset 200 is in 4 proteins but is not matched by any query.
        self.assertEqual(set(), gpi.offsets[200]['readIds'])
        expected = set(['AJF208%02d.1' % i for i in range(4, 8)])
        self.assertEqual(expected, gpi.offsets[200]['proteinAccessions'])

        # Offset 0 is in 3 proteins and is matched by query1.
        self.assertEqual({'query1'}, gpi.offsets[0]['readIds'])
        expected = set(['AJF208%02d.1' % i for i in range(4, 7)])
        self.assertEqual(expected, gpi.offsets[0]['proteinAccessions'])

        # Offset 1400 is in 2 proteins and is matched by query2.
        self.assertEqual({'query2'}, gpi.offsets[1400]['readIds'])
        self.assertEqual({'AJF20804.1', 'AJF20808.1'},
                         gpi.offsets[1400]['proteinAccessions'])

        # Read ids for all proteins.
        self.assertEqual({'query1', 'query2'}, gpi.readIdsForAllProteins())

        # AJF20804.1 coverage (its ranges are 2306-3221 and 0-1623)
        info = gpi.proteinCoverageInfo('AJF20804.1')
        self.assertEqual(300, info['coveredOffsets'])
        self.assertEqual(300, info['totalBases'])
        self.assertEqual((3221 - 2306) + (1623 - 0), info['ntLength'])
        self.assertEqual({'query1', 'query2'}, info['readIds'])

    def testLoadBAM123(self):
        """
        Test that everything is as expected after loading the BAM1, BAM2,
        and BAM3 files.
        """
        gpi = GenomeProteinInfo('KJ586809.1', DB, True)
        gpi.addSAM(BAM1)
        gpi.addSAM(BAM2)
        gpi.addSAM(BAM3)

        # Genome covered offsets.
        # There are 50 offsets that are covered twice.
        self.assertEqual(750 - 50, len(gpi.coveredOffsetCount))
        self.assertEqual(set(range(200)) | set(range(1000, 1500)),
                         set(gpi.coveredOffsetCount))
        self.assertEqual(set([1] * 700 + [2] * 50),
                         set(gpi.coveredOffsetCount.values()))

        # SAM files.
        self.assertEqual([BAM1, BAM2, BAM3], gpi.samFiles)
        self.assertEqual({'query1', 'query2', 'query3'},
                         gpi.readIdsMatchingGenome)

        # Protein accession numbers.
        expected = set(['AJF208%02d.1' % i for i in range(4, 11)])
        self.assertEqual(expected, set(gpi.proteins))

        # Offset 200 is in 4 proteins but is not matched by any query.
        self.assertEqual(set(), gpi.offsets[200]['readIds'])
        expected = set(['AJF208%02d.1' % i for i in range(4, 8)])
        self.assertEqual(expected, gpi.offsets[200]['proteinAccessions'])

        # Offset 0 is in 3 proteins and is matched by query1.
        self.assertEqual({'query1'}, gpi.offsets[0]['readIds'])
        expected = set(['AJF208%02d.1' % i for i in range(4, 7)])
        self.assertEqual(expected, gpi.offsets[0]['proteinAccessions'])

        # Offset 1400 is in 2 proteins and is matched by query2 and query3.
        self.assertEqual({'query2', 'query3'}, gpi.offsets[1400]['readIds'])
        self.assertEqual({'AJF20804.1', 'AJF20808.1'},
                         gpi.offsets[1400]['proteinAccessions'])

        # Read ids for all proteins.
        self.assertEqual({'query1', 'query2', 'query3'},
                         gpi.readIdsForAllProteins())

        # AJF20804.1 coverage (its ranges are 2306-3221 and 0-1623)
        info = gpi.proteinCoverageInfo('AJF20804.1')
        self.assertEqual(700, info['coveredOffsets'])
        self.assertEqual(750, info['totalBases'])
        self.assertEqual((3221 - 2306) + (1623 - 0), info['ntLength'])
        self.assertEqual({'query1', 'query2', 'query3'}, info['readIds'])

    def testTooFewReadOffsetsBAM1(self):
        """
        Test that a read is not returned as overlapping a protein unless it
        meets the minimum number of required overlapping offsets.
        """
        gpi = GenomeProteinInfo('KJ586809.1', DB, True)
        gpi.addSAM(BAM1)

        # Look at protein AJF20804.1 coverage (its ranges are 2306-3221 and
        # 0-1623). There should be no matching reads because the query
        # (query1) is only 200 nt long and so cannot match with at least
        # 500 nucleotides. The number of covered offsets and total bases
        # should both also be zero for the same reason.
        info = gpi.proteinCoverageInfo('AJF20804.1', 500)
        self.assertEqual(set(), info['readIds'])
        self.assertEqual(0, info['totalBases'])
        self.assertEqual(0, info['coveredOffsets'])

    def testSufficientReadOffsetsBAM1(self):
        """
        Test that a read is returned as overlapping a protein when it meets
        the minimum number of required overlapping offsets.
        """
        gpi = GenomeProteinInfo('KJ586809.1', DB, True)
        gpi.addSAM(BAM1)

        # Look at protein AJF20804.1 coverage (its ranges are 2306-3221 and
        # 0-1623). The query (query1) must be returned as it has 200
        # matching nucleotides.
        info = gpi.proteinCoverageInfo('AJF20804.1', 199)
        self.assertEqual({'query1'}, info['readIds'])

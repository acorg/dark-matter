from six import assertRaisesRegex
from unittest import TestCase
from tempfile import mkstemp
from os import close, unlink, write
from contextlib import contextmanager

from dark.reads import Read
from dark.sam import PaddedSAM, UnequalReferenceLengthError, UnknownReference


# These tests actually use the filesystem to read files. That's due to the API
# to pysam and the fact that it calls a C function to open files, so we can't
# mock Python's 'open' method. Hence the following context manager.
@contextmanager
def dataFile(data):
    """
    Create a context manager to store data in a temporary file and
    later remove it.
    """
    fd, filename = mkstemp()
    write(fd, data.encode('utf-8'))
    close(fd)
    yield filename
    unlink(filename)


class TestPaddedSAM(TestCase):
    """
    Test the PaddedSAM class.
    """

    # In reading the tests below, it is important to remember that the start
    # position (in the reference) of the match in SAM format is 1-based.  This
    # is the 4th field in the non-header SAM lines (i.e., those that don't
    # start with @). If you look at the code in ../dark/sam.py, pysam provides
    # a 'reference_start' attribute that is 0-based.

    def testReferencesToStr(self):
        """
        The referencesToStr method must return the expected string.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            self.assertEqual('id1 (length 90)\nid2 (length 91)',
                             ps.referencesToStr())
            ps.close()

    def testUnknownReferences(self):
        """
        Passing an unknown reference name to 'queries' must result in an
        UnknownReference exception.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            error = ("^Reference 'unknown' is not present in the "
                     "SAM/BAM file\\.$")
            queries = ps.queries(referenceName='unknown')
            assertRaisesRegex(self, UnknownReference, error, list, queries)
            ps.close()

    def testUnequalReferenceLengths(self):
        """
        Passing no reference name to 'queries' when the references have
        different lengths must result in an UnequalReferenceLengthError
        exception.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            error = ('^Your SAM/BAM file has 2 reference sequences, and their '
                     'lengths \(90, 91\) are not all identical\.$')
            queries = ps.queries()
            assertRaisesRegex(self, UnequalReferenceLengthError, error, list,
                              queries)
            ps.close()

    def testAllMMatch(self):
        """
        A simple all-'M' match must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 6M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testMixedMatch(self):
        """
        A match that is a mix of M, =, and X must result in the expected
        padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testHardClipLeft(self):
        """
        A simple all-'M' match with a hard clip left must result in the
        expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 10H6M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testHardClipRight(self):
        """
        A simple all-'M' match with a hard clip right must result in the
        expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 6M10H * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testRcNeeded(self):
        """
        A reverse-complimented match (flag = 16) when rcNeeded=True is passed
        must result in the expected (reverse complimented) padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 16 ref1 2 60 6M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(rcNeeded=True))
            self.assertEqual(Read('query1', '-CCTAGA---'), read)
            ps.close()

    def testRcSuffix(self):
        """
        A reverse-complimented sequence should have the rcSuffix string added
        to its id when an rcSuffix value is passed.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 16 ref1 2 60 6M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(rcSuffix='-rc'))
            self.assertEqual(Read('query1-rc', '-TCTAGG---'), read)
            ps.close()

    def testQuerySoftClipLeft(self):
        """
        A match with a soft-clipped region that does not extend to the left of
        the reference must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 4 60 2S4M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testQuerySoftClipReachesLeftEdge(self):
        """
        A match with a soft-clipped region that reaches to the left edge of the
        reference must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 5 60 4S2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'TCTAGG----'), read)
            ps.close()

    def testQuerySoftClipProtrudesLeft(self):
        """
        A match with a soft-clipped region that extends to the left of the
        reference must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 4S2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'AGG-------'), read)
            ps.close()

    def testKF414679SoftClipLeft(self):
        """
        Test for a case that wasn't working.
        """
        seq = ('GCCATGCAGTGGAACTCCACAGCATTCCACCAAGCTCTGC'
               'AGAATCCCAAAGTCAGGGGTTTGTATCTTCTTGCTGGTGGC')
        quality = ('ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ'
                   'ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ')
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 5 60 18S63M * 0 0 %s %s' % (seq, quality),
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', seq[14:]), read)
            ps.close()

    def testQuerySoftClipRight(self):
        """
        A match with a soft-clipped region that does not extend to the right of
        the reference must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 4 60 4M2S * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '---TCTAGG-'), read)
            ps.close()

    def testQuerySoftClipReachesRightEdge(self):
        """
        A match with a soft-clipped region that reaches to the right edge of
        the reference must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 5 60 2M4S * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '----TCTAGG'), read)
            ps.close()

    def testQuerySoftClipProtrudesRight(self):
        """
        A match with a soft-clipped region that extends to the right of
        the reference must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 6 60 2M4S * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-----TCTAG'), read)
            ps.close()

    def testQuerySoftClipProtrudesBothSides(self):
        """
        A match with a soft-clipped region that extends to both the left and
        right of the reference must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 4 60 5S5M5S * 0 0 TCTAGGCTGACTAAG ZZZZZZZZZZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'TAGGCTGACT'), read)
            ps.close()

    def testQueryHardClipAndSoftClipProtrudesBothSides(self):
        """
        A match with a soft-clipped region that extends to both the left and
        right of the reference must result in the expected padded sequence
        when hard clipping is also indicated by the CIGAR string.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 4 0 3H5S5M4S5H * 0 0 TCTAGGCTGACTAA ZZZZZZZZZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'TAGGCTGACT'), read)
            ps.close()

    def testReferenceInsertion(self):
        """
        An insertion into the reference must result in the expected padded
        sequence and the expected value in the referenceInsertions dictionary.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2M2I2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCGG-----'), read)
            self.assertEqual(
                {
                    3: {'T': 1},
                    4: {'A': 1},
                },
                ps.referenceInsertions)
            ps.close()

    def testReferenceDeletion(self):
        """
        An deletion of reference bases must result in the expected padded
        sequence (with gaps).
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2M2D4M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCNNTAGG-'), read)
            ps.close()

    def testReferenceDeletionAlternateChar(self):
        """
        An deletion of reference bases must result in the expected padded
        sequence (with gaps) when a queryInsertionChar is passed
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2M2D4M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(queryInsertionChar='?'))
            self.assertEqual(Read('query1', '-TC??TAGG-'), read)
            ps.close()

    def testReferenceSkip(self):
        """
        An skip of reference bases must result in the expected padded
        sequence (with gaps).
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2M2N4M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCNNTAGG-'), read)
            ps.close()

    def testReferenceSkipAlternateChar(self):
        """
        An skip of reference bases must result in the expected padded
        sequence (with gaps) when a queryInsertionChar is passed.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2M2N4M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(queryInsertionChar='X'))
            self.assertEqual(Read('query1', '-TCXXTAGG-'), read)
            ps.close()

    def testMixedMatchSpecificReferenceButNoMatches(self):
        """
        A request for reads aligned against a reference that exists but that
        has no matches must result in an empty list.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            '@SQ SN:ref2 LN:15',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            self.assertEqual([], list(ps.queries(referenceName='ref2')))
            ps.close()

    def testMixedMatchSpecificReference(self):
        """
        A match that is a mix of M, =, and X must result in the expected
        padded sequence when a reference sequence is specified.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            '@SQ SN:ref2 LN:15',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(referenceName='ref1'))
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testMinLength(self):
        """
        A request for reads longer than a certain value should result
        in the expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 0 ref1 2 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(minLength=6))
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testMinLengthWithReferenceDeletion(self):
        """
        The minLength specification must be applied after deletion of
        reference bases (which results in the query being lengthened to
        continue the match).
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2M2D4M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(minLength=7))
            self.assertEqual(Read('query1', '-TCNNTAGG-'), read)
            ps.close()

    def testDropSecondary(self):
        """
        Dropping matches flagged as secondary must give the expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 256 ref1 2 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(dropSecondary=True))
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testDropSupplementary(self):
        """
        Dropping matches flagged as supplementary must give the expected
        result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 2048 ref1 2 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(dropSupplementary=True))
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testDropDuplicates(self):
        """
        Dropping matches flagged as optical or PCR duplicates must give the
        expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 1024 ref1 2 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read,) = list(ps.queries(dropDuplicates=True))
            self.assertEqual(Read('query1', '-TCTAGG---'), read)
            ps.close()

    def testAllowDuplicateIds(self):
        """
        It must be possible to allow duplicate ids (in this case due to a
        secondary match).
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query1 0 ref1 3 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read1, read2) = list(ps.queries(allowDuplicateIds=True))
            self.assertEqual(Read('query1', '-TCTAGG---'), read1)
            self.assertEqual(Read('query1', '--TC------'), read2)
            ps.close()

    def testDuplicateIdDisambiguation(self):
        """
        Duplicate ids must be disambiguated if allowDuplicateIds is not given.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query1 0 ref1 3 60 2= * 0 0 TC ZZ',
            'query1 0 ref1 5 60 2S2= * 0 0 TCGA ZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read1, read2, read3) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read1)
            self.assertEqual(Read('query1/1', '--TC------'), read2)
            self.assertEqual(Read('query1/2', '--TCGA----'), read3)
            ps.close()

    def testKeepQualityControlFailures(self):
        """
        Keeping matches flagged as quality control failures must give the
        expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 512 ref1 4 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(filename)
            (read1, read2) = list(ps.queries(keepQCFailures=True))
            self.assertEqual(Read('query1', '-TCTAGG---'), read1)
            self.assertEqual(Read('query2', '---TC-----'), read2)
            ps.close()

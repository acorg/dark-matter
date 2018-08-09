from six import assertRaisesRegex
from unittest import TestCase
from tempfile import mkstemp
from os import close, unlink, write
from contextlib import contextmanager

from dark.reads import Read, ReadFilter
from dark.sam import (
    PaddedSAM, SAMFilter, UnequalReferenceLengthError, UnknownReference,
    InvalidSAM, samReferencesToStr)


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


class TestSAMFilter(TestCase):
    """
    Test the SAMFilter class.
    """
    def testUnknownReferences(self):
        """
        Passing an unknown reference id to the referenceLengths method must
        result in an UnknownReference exception.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:90',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sam = SAMFilter(filename, referenceIds={'unknown'})
            error = ("^Reference 'unknown' is not present in the "
                     "SAM/BAM file\\.$")
            assertRaisesRegex(self, UnknownReference, error,
                              sam.referenceLengths)

    def testStoreQueryIds(self):
        """
        If we request that query ids are saved, they must be.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456',
            'query2 0 ref1 2 60 2= * 0 0 TC XY',
            'query2 0 ref1 2 60 2= * 0 0 TC XY',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sf = SAMFilter(filename, storeQueryIds=True)
            list(sf.alignments())
            self.assertEqual({'query1', 'query2'}, sf.queryIds)

    def testAlignmentCount(self):
        """
        When all queries have been yielded, the alignment count must be
        as expected.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456',
            'query2 0 ref1 2 60 2= * 0 0 TC XY',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sf = SAMFilter(filename)
            list(sf.alignments())
            self.assertEqual(2, sf.alignmentCount)

    def testMinLength(self):
        """
        A request for reads that are only longer than a certain value should
        result in the expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 0 ref1 2 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            filterRead = ReadFilter(minLength=6).filter
            sf = SAMFilter(filename, filterRead=filterRead)
            (alignment,) = list(sf.alignments())
            self.assertEqual('query1', alignment.query_name)

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
            sf = SAMFilter(filename, dropSecondary=True)
            (alignment,) = list(sf.alignments())
            self.assertEqual('query1', alignment.query_name)

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
            sf = SAMFilter(filename, dropSupplementary=True)
            (alignment,) = list(sf.alignments())
            self.assertEqual('query1', alignment.query_name)

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
            sf = SAMFilter(filename, dropDuplicates=True)
            (alignment,) = list(sf.alignments())
            self.assertEqual('query1', alignment.query_name)

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
            sf = SAMFilter(filename, keepQCFailures=True)
            (alignment1, alignment2) = list(sf.alignments())
            self.assertEqual('query1', alignment1.query_name)
            self.assertEqual('query2', alignment2.query_name)

    def testMinScoreNoScores(self):
        """
        A request for reads with alignment scores no lower than a given value
        must produce an empty result when no alignments have scores.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 0 ref1 2 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sf = SAMFilter(filename, minScore=6)
            self.assertEqual([], list(sf.alignments()))

    def testMinScore(self):
        """
        A request for reads with alignment scores no lower than a given value
        must produce the expected result when some alignments have scores.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:10',
            'query2 0 ref1 2 60 2= * 0 0 TC ZZ',
            'query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:3',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sf = SAMFilter(filename, minScore=6)
            (alignment,) = list(sf.alignments())
            self.assertEqual('query1', alignment.query_name)

    def testMaxScoreNoScores(self):
        """
        A request for reads with alignment scores no higher than a given value
        must produce an empty result when no alignments have scores.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 0 ref1 2 60 2= * 0 0 TC ZZ',
            'query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sf = SAMFilter(filename, maxScore=6)
            self.assertEqual([], list(sf.alignments()))

    def testMaxScore(self):
        """
        A request for reads with alignment scores no higher than a given value
        must produce the expected result when some alignments have scores.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:10',
            'query2 0 ref1 2 60 2= * 0 0 TC ZZ',
            'query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:3',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sf = SAMFilter(filename, maxScore=6)
            (alignment,) = list(sf.alignments())
            self.assertEqual('query3', alignment.query_name)

    def testMinAndMaxScore(self):
        """
        A request for reads with alignment scores no lower or higher than
        given values must produce the expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:10',
            'query2 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:12',
            'query3 0 ref1 2 60 2= * 0 0 TC ZZ',
            'query4 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:3',
            'query5 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:2',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            sf = SAMFilter(filename, minScore=3, maxScore=10)
            (alignment1, alignment2) = list(sf.alignments())
            self.assertEqual('query1', alignment1.query_name)
            self.assertEqual('query4', alignment2.query_name)


class TestPaddedSAM(TestCase):
    """
    Test the PaddedSAM class.
    """

    # In reading the tests below, it is important to remember that the start
    # position (in the reference) of the match in SAM format is 1-based.  This
    # is the 4th field in the non-header SAM lines (i.e., those that don't
    # start with @). If you look at the code in ../dark/sam.py, pysam provides
    # a 'reference_start' attribute that is 0-based.

    def testUnequalReferenceLengths(self):
        """
        Passing no reference ids when the references have different lengths
        must result in an UnequalReferenceLengthError exception.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            error = ('^Your 2 SAM/BAM file reference sequence lengths '
                     '\\(id1=90, id2=91\\) are not all identical\\.$')
            assertRaisesRegex(self, UnequalReferenceLengthError, error,
                              PaddedSAM, SAMFilter(filename))

    def testAllMMatch(self):
        """
        A simple all-'M' match must result in the expected padded sequence.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 6M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(rcNeeded=True))
            self.assertEqual(Read('query1', '-CCTAGA---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(rcSuffix='-rc'))
            self.assertEqual(Read('query1-rc', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'TCTAGG----'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'AGG-------'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', seq[14:]), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '---TCTAGG-'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '----TCTAGG'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-----TCTAG'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'TAGGCTGACT'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', 'TAGGCTGACT'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCGG-----'), read)
            self.assertEqual(
                {
                    'query1': [(3, 'TA')],
                },
                ps.referenceInsertions)

    def testPrimaryAndSecondaryReferenceInsertion(self):
        """
        A primary and secondary insertion into the reference (of the same
        query) must result in the expected padded sequences and the expected
        value in the referenceInsertions dictionary.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2M2I2M * 0 0 TCTAGG ZZZZZZ',
            'query1 256 ref1 4 60 2M3I1M * 0 0 * *',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCGG-----'), read1)
            self.assertEqual(Read('query1/1', '---TCG----'), read2)
            self.assertEqual(
                {
                    'query1': [(3, 'TA')],
                    'query1/1': [(5, 'TAG')],
                },
                ps.referenceInsertions)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCNNTAGG-'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(queryInsertionChar='?'))
            self.assertEqual(Read('query1', '-TC??TAGG-'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCNNTAGG-'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(queryInsertionChar='X'))
            self.assertEqual(Read('query1', '-TCXXTAGG-'), read)

    def testMixedMatchSpecificReferenceButNoMatches(self):
        """
        A request for reads aligned against a reference that exists but that
        has no matches must result in an empty list.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:15',
            '@SQ SN:ref2 LN:15',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, referenceIds={'ref2'}))
            self.assertEqual([], list(ps.queries()))

    def testMixedMatchSpecificReference(self):
        """
        A match that is a mix of M, =, and X must result in the expected
        padded sequence when a reference sequence is specified.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            '@SQ SN:ref2 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, referenceIds={'ref1'}))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

    def testMinLength(self):
        """
        A request for reads that are only longer than a certain value should
        result in the expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ',
            'query2 0 ref1 2 60 2= * 0 0 TC ZZ',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            filterRead = ReadFilter(minLength=6).filter
            ps = PaddedSAM(SAMFilter(filename, filterRead=filterRead))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename, dropSecondary=True))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename, dropSupplementary=True))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename, dropDuplicates=True))
            (read,) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries(allowDuplicateIds=True))
            self.assertEqual(Read('query1', '-TCTAGG---'), read1)
            self.assertEqual(Read('query1', '--TC------'), read2)

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
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2, read3) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read1)
            self.assertEqual(Read('query1/1', '--TC------'), read2)
            self.assertEqual(Read('query1/2', '--TCGA----'), read3)

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
            ps = PaddedSAM(SAMFilter(filename, keepQCFailures=True))
            (read1, read2) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCTAGG---'), read1)
            self.assertEqual(Read('query2', '---TC-----'), read2)

    def testSecondaryWithNoPreviousSequence(self):
        """
        A secondary match with a '*' seq that is not preceeded by a query with
        a sequence must result in a ValueError being raised.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query 256 ref1 3 60 4M * 0 0 * *',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            error = ('^Query line 1 has an empty SEQ field, but no previous '
                     'alignment is present\\.$')
            queries = ps.queries()
            assertRaisesRegex(self, InvalidSAM, error, list, queries)

    def testSecondaryWithNoSequence(self):
        """
        A secondary match with a '*' seq must result in the sequence from the
        previous query being used.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 3M * 0 0 TCT ZZZ',
            'query2 0 ref1 2 60 4M * 0 0 TCTA ZZZZ',
            'query2 256 ref1 6 60 4M * 0 0 * *',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2, read3) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCT------'), read1)
            self.assertEqual(Read('query2', '-TCTA-----'), read2)
            self.assertEqual(Read('query2/1', '-----TCTA-'), read3)

    def testSupplementaryWithNoPreviousSequence(self):
        """
        A supplementary match with a '*' seq that is not preceeded by a query
        with a sequence must result in a ValueError being raised.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query 2048 ref1 3 60 4M * 0 0 * *',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            error = ('^Query line 1 has an empty SEQ field, but no previous '
                     'alignment is present\\.$')
            queries = ps.queries()
            assertRaisesRegex(self, InvalidSAM, error, list, queries)

    def testSupplementaryWithNoSequence(self):
        """
        A supplementary match with a '*' seq must result in the sequence from
        the previous query being used.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 3M * 0 0 TCT ZZZ',
            'query2 0 ref1 2 60 4M * 0 0 TCTA ZZZZ',
            'query2 2048 ref1 6 60 4M * 0 0 * *',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2, read3) = list(ps.queries())
            self.assertEqual(Read('query1', '-TCT------'), read1)
            self.assertEqual(Read('query2', '-TCTA-----'), read2)
            self.assertEqual(Read('query2/1', '-----TCTA-'), read3)

    def testNotSecondaryAndNotSupplementaryWithNoSequence(self):
        """
        An alignment with a '*' seq that is not secondary or supplementary
        must result in a ValueError being raised.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query 0 ref1 3 60 4M * 0 0 * *',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            error = ('^Query line 1 has an empty SEQ field, but the alignment '
                     'is not marked as secondary or supplementary\\.$')
            queries = ps.queries()
            assertRaisesRegex(self, InvalidSAM, error, list, queries)

    def testAlsoYieldAlignments(self):
        """
        A request for queries with their pysam alignments should have the
        expected result.
        """
        data = '\n'.join([
            '@SQ SN:ref1 LN:10',
            'query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456',
            'query2 0 ref1 2 60 2= * 0 0 TC XY',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries(addAlignment=True))

            self.assertEqual(Read('query1', '-TCTAGG---'), read1)
            self.assertEqual('TCTAGG', read1.alignment.query_sequence)
            self.assertEqual('123456', ''.join(
                map(lambda x: chr(x + 33), read1.alignment.query_qualities)))

            self.assertEqual(Read('query2', '-TC-------'), read2)
            self.assertEqual('TC', read2.alignment.query_sequence)
            self.assertEqual('XY', ''.join(
                map(lambda x: chr(x + 33), read2.alignment.query_qualities)))


class TestSamReferencesToStr(TestCase):
    """
    Test the samReferencesToStr function.
    """
    def testSimple(self):
        """
        The referencesToStr method must return the expected string.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            self.assertEqual('id1 (length 90)\nid2 (length 91)',
                             samReferencesToStr(filename))

    def testIndent(self):
        """
        The referencesToStr method must return the expected string when
        passed an indent.
        """
        data = '\n'.join([
            '@SQ SN:id1 LN:90',
            '@SQ SN:id2 LN:91',
        ]).replace(' ', '\t')

        with dataFile(data) as filename:
            self.assertEqual('  id1 (length 90)\n  id2 (length 91)',
                             samReferencesToStr(filename, indent=2))

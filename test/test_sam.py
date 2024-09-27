from unittest import TestCase
from tempfile import mkstemp
from os import close, unlink, write
from contextlib import contextmanager
from io import StringIO
from json import dumps, loads
import numpy as np
from dataclasses import dataclass

from pysam import CHARD_CLIP, CMATCH

from dark.reads import Read, ReadFilter
from dark.sam import (
    DistanceMatrix,
    InvalidSAM,
    PaddedSAM,
    ReferenceReads,
    SAMFilter,
    UnequalReferenceLengthError,
    UnknownReference,
    _hardClip,
    samReferencesToStr,
)

# Note: getReferenceInfo in dark/sam.py is tested by the calls to
# consensusFromBAM in test/test_consensus.py


# Some tests below actually use the filesystem to read files. That's due to
# the API to pysam and the fact that it calls a C function to open files,
# so we can't mock Python's 'open' method. Hence the following context
# manager.
@contextmanager
def dataFile(data):
    """
    Create a context manager to store data in a temporary file and
    later remove it.
    """
    fd, filename = mkstemp()
    write(fd, data.encode("utf-8"))
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
        data = "\n".join(
            [
                "@SQ SN:id1 LN:90",
                "@SQ SN:id2 LN:90",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sam = SAMFilter(filename, referenceIds={"unknown"})
            error = "^Reference 'unknown' is not present in the SAM/BAM file\\.$"
            self.assertRaisesRegex(UnknownReference, error, sam.referenceLengths)

    def testNoFilteringOptions(self):
        """
        If no filtering options are given, the noFiltering attribute
        on the SAM filter must be True.
        """
        sf = SAMFilter(None)
        self.assertTrue(sf.noFiltering)

    def testNoFilteringAllAlignmentsReturned(self):
        """
        When no filtering options are given, all alignments must be returned.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:10",
                "query2 0 ref1 2 60 2= * 0 0 TC ZZ",
                "query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:3",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename)
            (alignment1, alignment2, alignment3) = list(sf.alignments())
            self.assertEqual("query1", alignment1.query_name)
            self.assertEqual("query2", alignment2.query_name)
            self.assertEqual("query3", alignment3.query_name)

    def testAFilteringOptionSetsNoFiltering(self):
        """
        If a filtering options is given, the noFiltering attribute
        on the SAM filter must be False.
        """
        sf = SAMFilter(None, storeQueryIds=True)
        self.assertFalse(sf.noFiltering)

    def testStoreQueryIds(self):
        """
        If we request that query ids are saved, they must be.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref1 2 60 2= * 0 0 TC XY",
                "query2 0 ref1 2 60 2= * 0 0 TC XY",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, storeQueryIds=True)
            list(sf.alignments())
            self.assertEqual({"query1", "query2"}, sf.queryIds)

    def testAlignmentCount(self):
        """
        When all queries have been yielded, the alignment count must be
        as expected.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref1 2 60 2= * 0 0 TC XY",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename)
            list(sf.alignments())
            self.assertEqual(2, sf.alignmentCount)

    def testMinLength(self):
        """
        A request for reads that are only longer than a certain value should
        result in the expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 0 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            filterRead = ReadFilter(minLength=6).filter
            sf = SAMFilter(filename, filterRead=filterRead)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query1", alignment.query_name)

    def testDropSecondary(self):
        """
        Dropping matches flagged as secondary must give the expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 256 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, dropSecondary=True)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query1", alignment.query_name)

    def testDropSupplementary(self):
        """
        Dropping matches flagged as supplementary must give the expected
        result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 2048 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, dropSupplementary=True)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query1", alignment.query_name)

    def testDropDuplicates(self):
        """
        Dropping matches flagged as optical or PCR duplicates must give the
        expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 1024 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, dropDuplicates=True)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query1", alignment.query_name)

    def testKeepQualityControlFailures(self):
        """
        Keeping matches flagged as quality control failures must give the
        expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 512 ref1 4 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, keepQCFailures=True)
            (alignment1, alignment2) = list(sf.alignments())
            self.assertEqual("query1", alignment1.query_name)
            self.assertEqual("query2", alignment2.query_name)

    def testMinScoreNoScores(self):
        """
        A request for reads with alignment scores no lower than a given value
        must produce an empty result when no alignments have scores.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 0 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, minScore=6)
            self.assertEqual([], list(sf.alignments()))

    def testMinScore(self):
        """
        A request for reads with alignment scores no lower than a given value
        must produce the expected result when some alignments have scores.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:10",
                "query2 0 ref1 2 60 2= * 0 0 TC ZZ",
                "query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:3",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, minScore=6)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query1", alignment.query_name)

    def testMaxScoreNoScores(self):
        """
        A request for reads with alignment scores no higher than a given value
        must produce an empty result when no alignments have scores.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 0 ref1 2 60 2= * 0 0 TC ZZ",
                "query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, maxScore=6)
            self.assertEqual([], list(sf.alignments()))

    def testMaxScore(self):
        """
        A request for reads with alignment scores no higher than a given value
        must produce the expected result when some alignments have scores.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:10",
                "query2 0 ref1 2 60 2= * 0 0 TC ZZ",
                "query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:3",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, maxScore=6)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query3", alignment.query_name)

    def testMinAndMaxScore(self):
        """
        A request for reads with alignment scores no lower or higher than
        given values must produce the expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:10",
                "query2 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:12",
                "query3 0 ref1 2 60 2= * 0 0 TC ZZ",
                "query4 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:3",
                "query5 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ AS:i:2",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename, minScore=3, maxScore=10)
            (alignment1, alignment2) = list(sf.alignments())
            self.assertEqual("query1", alignment1.query_name)
            self.assertEqual("query4", alignment2.query_name)

    def testCloseButNoCIGAR(self):
        """
        An unmapped query with no CIGAR string must be passed through
        unchanged if dropUnmapped is not specified.
        """
        data = "\n".join(
            [
                "@SQ SN:ref LN:10",
                "query1 4 * 0 0 * * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query1", alignment.query_name)
            self.assertEqual("TCTAGG", alignment.query_sequence)
            self.assertEqual(
                "ZZZZZZ", "".join(map(lambda x: chr(x + 33), alignment.query_qualities))
            )

    def testNoQuality(self):
        """
        If an alignment has * for the quality string, the filter must
        return an alignment with a C{None} quality value.
        """
        data = "\n".join(
            [
                "@SQ SN:ref LN:10",
                "query1 4 * 0 0 6M * 0 0 TCTAGG *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            sf = SAMFilter(filename)
            (alignment,) = list(sf.alignments())
            self.assertEqual("query1", alignment.query_name)
            self.assertEqual("TCTAGG", alignment.query_sequence)
            self.assertIsNone(alignment.query_qualities)


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
        data = "\n".join(
            [
                "@SQ SN:id1 LN:90",
                "@SQ SN:id2 LN:91",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            error = (
                "^Your 2 SAM/BAM file reference sequence lengths "
                "\\(id1=90, id2=91\\) are not all identical\\.$"
            )
            self.assertRaisesRegex(
                UnequalReferenceLengthError, error, PaddedSAM, SAMFilter(filename)
            )

    def testQueryTooLong(self):
        """
        If the query sequence is longer than the total of the lengths in the
        CIGAR operations, a ValueError must be raised.
        """
        # This test just returns. It used to be possible to reach the
        # "Query ... not fully consumed when parsing CIGAR string."
        # ValueError in sam.py, prior to the fix of
        # https://github.com/acorg/dark-matter/issues/630 but it is not
        # possible to get a CIGAR string that has a different total length
        # from the sequence length through to our code in sam.py because
        # pysam catches the error.  I'm leaving this test here because it
        # documents that the error checked for in sam.py cannot currently
        # be reached and the test may become useful. For now it just returns.
        return
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:90",
                "query1 0 ref1 1 60 4M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            error = "^Query TCTAGG not fully consumed when parsing CIGAR " "string\\."
            self.assertRaisesRegex(ValueError, error, list, ps.queries())

    def testAllMMatch(self):
        """
        A simple all-'M' match must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 6M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testMixedMatch(self):
        """
        A match that is a mix of M, =, and X must result in the expected
        padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testHardClipLeft(self):
        """
        A simple all-'M' match with a hard clip left must result in the
        expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 10H6M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testHardClipRight(self):
        """
        A simple all-'M' match with a hard clip right must result in the
        expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 6M10H * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testRcNeeded(self):
        """
        A reverse-complemented match (flag = 16) when rcNeeded=True is passed
        must result in the expected (reverse complemented) padded sequence
        and reversed quality string.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 16 ref1 2 60 6M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(rcNeeded=True))
            self.assertEqual(Read("query1", "-CCTAGA---", "!654321!!!"), read)

    def testRcSuffix(self):
        """
        A reverse-complemented sequence should have the rcSuffix string added
        to its id when an rcSuffix value is passed.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 16 ref1 2 60 6M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(rcSuffix="-rc", rcNeeded=True))
            self.assertEqual(Read("query1-rc", "-CCTAGA---", "!654321!!!"), read)

    def testQuerySoftClipLeft(self):
        """
        A match with a soft-clipped region that does not extend to the left of
        the reference must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 4 60 2S4M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testQuerySoftClipReachesLeftEdge(self):
        """
        A match with a soft-clipped region that reaches to the left edge of the
        reference must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 5 60 4S2M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "TCTAGG----", "ZZZZZZ!!!!"), read)

    def testQuerySoftClipProtrudesLeft(self):
        """
        A match with a soft-clipped region that extends to the left of the
        reference must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 4S2M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "AGG-------", "ZZZ!!!!!!!"), read)

    def testKF414679SoftClipLeft(self):
        """
        Test for a case that wasn't working.
        """
        seq = (
            "GCCATGCAGTGGAACTCCACAGCATTCCACCAAGCTCTGC"
            "AGAATCCCAAAGTCAGGGGTTTGTATCTTCTTGCTGGTGGC"
        )
        quality = (
            "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
            "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
        )
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 5 60 18S63M * 0 0 %s %s" % (seq, quality),
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", seq[14:], quality[14:]), read)

    def testQuerySoftClipRight(self):
        """
        A match with a soft-clipped region that does not extend to the right of
        the reference must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 4 60 4M2S * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "---TCTAGG-", "!!!ZZZZZZ!"), read)

    def testQuerySoftClipReachesRightEdge(self):
        """
        A match with a soft-clipped region that reaches to the right edge of
        the reference must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 5 60 2M4S * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "----TCTAGG", "!!!!ZZZZZZ"), read)

    def testQuerySoftClipProtrudesRight(self):
        """
        A match with a soft-clipped region that extends to the right of
        the reference must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 6 60 2M4S * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-----TCTAG", "!!!!!ZZZZZ"), read)

    def testQuerySoftClipProtrudesBothSides(self):
        """
        A match with a soft-clipped region that extends to both the left and
        right of the reference must result in the expected padded sequence.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 4 60 5S5M5S * 0 0 TCTAGGCTGACTAAG ZZZZZZZZZZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "TAGGCTGACT", "ZZZZZZZZZZ"), read)

    def testQueryHardClipAndSoftClipProtrudesBothSides(self):
        """
        A match with a soft-clipped region that extends to both the left and
        right of the reference must result in the expected padded sequence
        when hard clipping is also indicated by the CIGAR string.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 4 0 3H5S5M4S5H * 0 0 TCTAGGCTGACTAA ZZZZZZZZZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "TAGGCTGACT", "ZZZZZZZZZZ"), read)

    def testReferenceInsertion(self):
        """
        An insertion into the reference must result in the expected padded
        sequence and the expected value in the referenceInsertions dictionary.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2M2I2M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCGG-----", "!ZZZZ!!!!!"), read)
            self.assertEqual(
                {
                    "query1": [(3, "TA")],
                },
                ps.referenceInsertions,
            )

    def testPrimaryAndSecondaryReferenceInsertion(self):
        """
        A primary and secondary insertion into the reference (of the same
        query) must result in the expected padded sequences and the expected
        value in the referenceInsertions dictionary.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2M2I2M * 0 0 TCTAGG ZZZZZZ",
                "query1 256 ref1 4 60 2M3I1M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCGG-----", "!ZZZZ!!!!!"), read1)
            self.assertEqual(Read("query1/1", "---TCG----", "!!!ZZZ!!!!"), read2)
            self.assertEqual(
                {
                    "query1": [(3, "TA")],
                    "query1/1": [(5, "TAG")],
                },
                ps.referenceInsertions,
            )

    def testReferenceDeletion(self):
        """
        An deletion of reference bases must result in the expected padded
        sequence (with Ns inserted for the deleted reference bases).
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2M2D4M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCNNTAGG-", "!ZZ!!ZZZZ!"), read)

    def testReferenceDeletionAlternateChars(self):
        """
        An deletion of reference bases must result in the expected padded
        sequence (with the passed query insertion character and unknown
        quality character) when queryInsertionChar and unknownQualityChar
        arguments are passed.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2M2D4M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(queryInsertionChar="?", unknownQualityChar="+"))
            self.assertEqual(Read("query1", "-TC??TAGG-", "+ZZ++ZZZZ+"), read)

    def testReferenceSkip(self):
        """
        An skip of reference bases must result in the expected padded
        sequence with the passed unknown quality character when the
        unknownQualityChar argument is passed.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2M2N4M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(unknownQualityChar="."))
            self.assertEqual(Read("query1", "-TCNNTAGG-", ".ZZ..ZZZZ."), read)

    def testReferenceSkipAlternateChars(self):
        """
        An skip of reference bases must result in the expected padded
        sequence (with the passed query insertion character and unknown
        quality character) when queryInsertionChar and unknownQualityChar
        arguments are passed.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2M2N4M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read,) = list(ps.queries(queryInsertionChar="X", unknownQualityChar="+"))
            self.assertEqual(Read("query1", "-TCXXTAGG-", "+ZZ++ZZZZ+"), read)

    def testMixedMatchSpecificReferenceButNoMatches(self):
        """
        A request for reads aligned against a reference that exists but that
        has no matches must result in an empty list.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:15",
                "@SQ SN:ref2 LN:15",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, referenceIds={"ref2"}))
            self.assertEqual([], list(ps.queries()))

    def testMixedMatchSpecificReference(self):
        """
        A match that is a mix of M, =, and X must result in the expected
        padded sequence when a reference sequence is specified.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, referenceIds={"ref1"}))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testMinLength(self):
        """
        A request for reads that are only longer than a certain value should
        result in the expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 0 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            filterRead = ReadFilter(minLength=6).filter
            ps = PaddedSAM(SAMFilter(filename, filterRead=filterRead))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testDropSecondary(self):
        """
        Dropping matches flagged as secondary must give the expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 256 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, dropSecondary=True))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testDropSupplementary(self):
        """
        Dropping matches flagged as supplementary must give the expected
        result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 2048 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, dropSupplementary=True))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testDropDuplicates(self):
        """
        Dropping matches flagged as optical or PCR duplicates must give the
        expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 1024 ref1 2 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, dropDuplicates=True))
            (read,) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read)

    def testAllowDuplicateIds(self):
        """
        It must be possible to allow duplicate ids (in this case due to a
        secondary match).
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query1 0 ref1 3 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries(allowDuplicateIds=True))
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read1)
            self.assertEqual(Read("query1", "--TC------", "!!ZZ!!!!!!"), read2)

    def testDuplicateIdDisambiguation(self):
        """
        Duplicate ids must be disambiguated if allowDuplicateIds is not given.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query1 0 ref1 3 60 2= * 0 0 TC ZZ",
                "query1 0 ref1 5 60 2S2= * 0 0 TCGA ZZZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2, read3) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read1)
            self.assertEqual(Read("query1/1", "--TC------", "!!ZZ!!!!!!"), read2)
            self.assertEqual(Read("query1/2", "--TCGA----", "!!ZZZZ!!!!"), read3)

    def testKeepQualityControlFailures(self):
        """
        Keeping matches flagged as quality control failures must give the
        expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query2 512 ref1 4 60 2= * 0 0 TC ZZ",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename, keepQCFailures=True))
            (read1, read2) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read1)
            self.assertEqual(Read("query2", "---TC-----", "!!!ZZ!!!!!"), read2)

    def testSecondaryWithNoPreviousSequence(self):
        """
        A secondary match with a '*' seq that is not preceeded by a query with
        a sequence must result in a ValueError being raised.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query 256 ref1 3 60 4M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            error = (
                "^pysam produced an alignment \\(number 1\\) with no "
                "query sequence without previously giving an alignment "
                "with a sequence\\.$"
            )
            queries = ps.queries()
            self.assertRaisesRegex(InvalidSAM, error, list, queries)

    def testSecondaryWithNoSequence(self):
        """
        A secondary match with a '*' seq must result in the sequence from the
        previous query being used.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 3M * 0 0 TCT ZZZ",
                "query2 0 ref1 2 60 4M * 0 0 TCTA ZZZZ",
                "query2 256 ref1 6 60 4M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2, read3) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCT------", "!ZZZ!!!!!!"), read1)
            self.assertEqual(Read("query2", "-TCTA-----", "!ZZZZ!!!!!"), read2)
            self.assertEqual(Read("query2/1", "-----TCTA-", "!!!!!ZZZZ!"), read3)

    def testSupplementaryWithNoPreviousSequence(self):
        """
        A supplementary match with a '*' seq that is not preceeded by a query
        with a sequence must result in a ValueError being raised.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query 2048 ref1 3 60 4M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            error = (
                "^pysam produced an alignment \\(number 1\\) with no "
                "query sequence without previously giving an alignment "
                "with a sequence\\.$"
            )
            queries = ps.queries()
            self.assertRaisesRegex(InvalidSAM, error, list, queries)

    def testSupplementaryWithNoSequence(self):
        """
        A supplementary match with a '*' seq must result in the sequence from
        the previous query being used.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 3M * 0 0 TCT ZZZ",
                "query2 0 ref1 2 60 4M * 0 0 TCTA ZZZZ",
                "query2 2048 ref1 6 60 4M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2, read3) = list(ps.queries())
            self.assertEqual(Read("query1", "-TCT------", "!ZZZ!!!!!!"), read1)
            self.assertEqual(Read("query2", "-TCTA-----", "!ZZZZ!!!!!"), read2)
            self.assertEqual(Read("query2/1", "-----TCTA-", "!!!!!ZZZZ!"), read3)

    def testNotSecondaryAndNotSupplementaryWithNoSequence(self):
        """
        An alignment with a '*' seq that is not secondary or supplementary
        must result in a ValueError being raised.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query 0 ref1 3 60 4M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            error = (
                "^pysam produced an alignment \\(number 1\\) with no "
                "query sequence without previously giving an alignment "
                "with a sequence\\.$"
            )
            queries = ps.queries()
            self.assertRaisesRegex(InvalidSAM, error, list, queries)

    def testAlsoYieldAlignments(self):
        """
        A request for queries with their pysam alignments should have the
        expected result.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref1 2 60 2= * 0 0 TC 78",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries(addAlignment=True))

            self.assertEqual(Read("query1", "-TCTAGG---", "!123456!!!"), read1)
            self.assertEqual("TCTAGG", read1.alignment.query_sequence)
            self.assertEqual(
                "123456",
                "".join(map(lambda x: chr(x + 33), read1.alignment.query_qualities)),
            )

            self.assertEqual(Read("query2", "-TC-------", "!78!!!!!!!"), read2)
            self.assertEqual("TC", read2.alignment.query_sequence)
            self.assertEqual(
                "78",
                "".join(map(lambda x: chr(x + 33), read2.alignment.query_qualities)),
            )

    def testHardClippingInCIGARButQueryNotHardClipped(self):
        """
        As documented in https://github.com/acorg/dark-matter/issues/630 we
        must deal correctly with a case in which the CIGAR string says a
        query is hard-clipped but the query sequence in the SAM file
        actually isn't. This can be due to a prior alignment with a soft clip,
        in which case the full query sequence has to be given before the
        secondary alignment with the hard clip.
        """
        data = "\n".join(
            [
                "@SQ SN:Chimp-D00220 LN:8",
                "@SQ SN:D-AM494716 LN:8",
                "@SQ SN:D-XXX LN:8",
                "@SQ SN:Chimp-YYY LN:8",
                "query1 0 Chimp-D00220 1 0 3S5M * 0 0 TTTTGGTT 12345678",
                "query1 256 D-AM494716 1 0 3H5M * 0 0 * *",
                "query1 256 D-XXX 1 0 5H3M * 0 0 * *",
                "query1 0 Chimp-YYY 1 0 8M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2, read3, read4) = list(ps.queries(addAlignment=True))

            self.assertEqual(Read("query1", "TGGTT---", "45678!!!"), read1)
            self.assertEqual("TTTTGGTT", read1.alignment.query_sequence)

            self.assertEqual(Read("query1/1", "TGGTT---", "45678!!!"), read2)
            self.assertEqual("TGGTT", read2.alignment.query_sequence)

            self.assertEqual(Read("query1/2", "GTT-----", "678!!!!!"), read3)
            self.assertEqual("GTT", read3.alignment.query_sequence)

            self.assertEqual(Read("query1/3", "TTTTGGTT", "12345678"), read4)
            self.assertEqual("TTTTGGTT", read4.alignment.query_sequence)

    def testSecondaryAlignmentHasQuery(self):
        """
        If the first alignment of a query is against a reference that is not
        wanted, a subsequent secondary alignment (SAM flag = 256) must have
        the original query and quality strings (even though these are only
        present in the SAM as * characters and the query is None when it comes
        back from pysam).
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query1 256 ref2 2 60 2=2X2M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries(addAlignment=True))
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read1)
            self.assertEqual("ref1", read1.alignment.reference_name)
            self.assertEqual(Read("query1/1", "-TCTAGG---", "!ZZZZZZ!!!"), read2)
            self.assertEqual("ref2", read2.alignment.reference_name)

    def testSupplementaryAlignmentHasQuery(self):
        """
        If the first alignment of a query is against a reference that is not
        wanted, a subsequent supplementary alignment (SAM flag = 2048) must
        have the original query and quality strings (even though these are only
        present in the SAM as * characters and the query is None when it comes
        back from pysam).
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG ZZZZZZ",
                "query1 2048 ref2 2 60 2=2X2M * 0 0 * *",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            ps = PaddedSAM(SAMFilter(filename))
            (read1, read2) = list(ps.queries(addAlignment=True))
            self.assertEqual(Read("query1", "-TCTAGG---", "!ZZZZZZ!!!"), read1)
            self.assertEqual("ref1", read1.alignment.reference_name)
            self.assertEqual(Read("query1/1", "-TCTAGG---", "!ZZZZZZ!!!"), read2)
            self.assertEqual("ref2", read2.alignment.reference_name)


class TestSamReferencesToStr(TestCase):
    """
    Test the samReferencesToStr function.
    """

    def testSimple(self):
        """
        The referencesToStr method must return the expected string.
        """
        data = "\n".join(
            [
                "@SQ SN:id1 LN:90",
                "@SQ SN:id2 LN:91",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            self.assertEqual(
                "id1 (length 90)\nid2 (length 91)", samReferencesToStr(filename)
            )

    def testIndent(self):
        """
        The referencesToStr method must return the expected string when
        passed an indent.
        """
        data = "\n".join(
            [
                "@SQ SN:id1 LN:90",
                "@SQ SN:id2 LN:91",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            self.assertEqual(
                "  id1 (length 90)\n  id2 (length 91)",
                samReferencesToStr(filename, indent=2),
            )


class TestHardClip(TestCase):
    """
    Test the _hardClip function.
    """

    def testHardClipInMiddle(self):
        """
        If hard clipping is given as an operation not at the beginning or end
        of the sequence, a ValueError must be raised.
        """
        error = (
            "^Invalid CIGAR tuples .* contains hard-clipping operation "
            "that is neither at the start nor the end of the sequence\\.$"
        )
        self.assertRaisesRegex(
            ValueError,
            error,
            _hardClip,
            "CGT",
            "123",
            (
                (CMATCH, 1),
                (CHARD_CLIP, 1),
                (CMATCH, 1),
            ),
        )

    def testThreeHardClips(self):
        """
        If hard clipping is specified more than twice, a ValueError must be
        raised.
        """
        error = (
            "^Invalid CIGAR tuples .* specifies hard-clipping 3 times "
            "\\(2 is the maximum\\).$"
        )
        self.assertRaisesRegex(
            ValueError,
            error,
            _hardClip,
            "CGT",
            "123",
            (
                (CHARD_CLIP, 1),
                (CHARD_CLIP, 1),
                (CHARD_CLIP, 1),
            ),
        )

    def testNoClip(self):
        """
        If no hard clipping is indicated, the function must return the
        original sequence.
        """
        self.assertEqual(("CGT", "123", False), _hardClip("CGT", "123", ((CMATCH, 3),)))

    def testClipLeft(self):
        """
        If hard clipping on the left is indicated, and has not been done,
        the function must return the expected sequence.
        """
        self.assertEqual(
            ("CGT", "456", True),
            _hardClip(
                "CAACGT",
                "123456",
                (
                    (CHARD_CLIP, 3),
                    (CMATCH, 3),
                ),
            ),
        )

    def testClipRight(self):
        """
        If hard clipping on the right is indicated, and has not been done,
        the function must return the expected sequence.
        """
        self.assertEqual(
            ("CA", "12", True),
            _hardClip(
                "CAACGT",
                "123456",
                (
                    (CMATCH, 2),
                    (CHARD_CLIP, 4),
                ),
            ),
        )

    def testClipBoth(self):
        """
        If hard clipping on the left and right is indicated, and has not been
        done, the function must return the expected sequence.
        """
        self.assertEqual(
            ("AA", "23", True),
            _hardClip(
                "CAACGT",
                "123456",
                (
                    (CHARD_CLIP, 1),
                    (CMATCH, 2),
                    (CHARD_CLIP, 3),
                ),
            ),
        )

    def testClipLeftAlreadyDone(self):
        """
        If hard clipping on the left is indicated, and has already been done,
        the function must return the expected sequence.
        """
        self.assertEqual(
            ("CGT", "123", False),
            _hardClip(
                "CGT",
                "123",
                (
                    (CHARD_CLIP, 3),
                    (CMATCH, 3),
                ),
            ),
        )

    def testClipRightAlreadyDone(self):
        """
        If hard clipping on the right is indicated, and has already been done,
        the function must return the expected sequence.
        """
        self.assertEqual(
            ("CA", "12", False),
            _hardClip(
                "CA",
                "12",
                (
                    (CMATCH, 2),
                    (CHARD_CLIP, 4),
                ),
            ),
        )

    def testClipBothAlreadyDone(self):
        """
        If hard clipping on the left and right is indicated, and has already
        been done, the function must return the expected sequence.
        """
        self.assertEqual(
            ("AA", "12", False),
            _hardClip(
                "AA",
                "12",
                (
                    (CHARD_CLIP, 1),
                    (CMATCH, 2),
                    (CHARD_CLIP, 3),
                ),
            ),
        )


class TestDistanceMatrix(TestCase):
    """
    Test the DistanceMatrix class.
    """

    def testEmpty(self):
        """
        The similarity (scores) matrix must be empty after intitialization.
        """
        dm = DistanceMatrix()
        self.assertEqual({}, dm.scores)

    def testNoMatchesHasNoMatchedReferenceIds(self):
        """
        If a SAM file with no query records is added, there should be no
        matched reference ids (even though there is a reference in the SAM
        header, it has no matches).
        """
        dm = DistanceMatrix()
        data = "\n".join(
            [
                "@SQ SN:id1 LN:90",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual({}, dm.scores)

    def testAddEmptySAMNoScoreTag(self):
        """
        If a SAM file with no query records is added and no score tag is
        passed, the scores matrix must be empty and the distance between
        two (non-existent) references must be 1.0.
        """
        dm = DistanceMatrix()
        data = "\n".join(
            [
                "@SQ SN:id1 LN:90",
            ]
        ).replace(" ", "\t")

        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual({}, dm.scores)
        self.assertEqual(1.0, dm.jaccardDistance("ref1", "ref2"))
        self.assertEqual(1.0, dm.soergelDistance("ref1", "ref2"))

    def testAddEmptySAMWithScoreTag(self):
        """
        If a SAM file with no query records is added and a score tag is passed,
        the similarity (scores) matrix must be empty.
        """
        data = "\n".join(
            [
                "@SQ SN:id1 LN:90",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual({}, dm.scores)
        self.assertEqual(1.0, dm.jaccardDistance("ref1", "ref2"))
        self.assertEqual(1.0, dm.soergelDistance("ref1", "ref2"))

    def testOneQueryMappedNoScoreTag(self):
        """
        If one query is mapped to one reference, the scores matrix must have
        a 1.0 score if no score tag is passed.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual(
            {
                "ref1": {
                    "query1": 1.0,
                },
            },
            dm.scores,
        )

    def testOneQueryMappedWithScoreTag(self):
        """
        If one query is mapped to one reference, the scores matrix must have
        the correct score if a score tag is passed.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:77",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(
            {
                "ref1": {
                    "query1": 77.0,
                },
            },
            dm.scores,
        )

        self.assertEqual(77.0, dm.score("ref1", "query1"))

    def testOneQueryMappedWithScoreTagFloat(self):
        """
        If one query is mapped to one reference, the scores matrix must have
        the correct score if a score tag is passed and the score is of type
        float (the AS:f:77.5 in the SAM record).
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:f:77.5",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(
            {
                "ref1": {
                    "query1": 77.5,
                },
            },
            dm.scores,
        )

        self.assertEqual(77.5, dm.score("ref1", "query1"))

    def testNonExistentQueryNotMapped(self):
        """
        If a query (not even existing in this case) is not mapped to the
        reference, the score between the two must be zero.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(0.0, dm.score("ref1", "query1"))

    def testNonExistentQuery(self):
        """
        The score for a non-existent query must be zero.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(0.0, dm.score("ref1", "query1"))

    def testQueryNotMapped(self):
        """
        If a query did not map to a reference, the score between the two must
        be zero.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:77",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(
            {
                "ref2": {
                    "query1": 77,
                },
            },
            dm.scores,
        )

        self.assertEqual(0.0, dm.score("ref1", "query1"))

    def testJaccardDistanceToSelf(self):
        """
        The Jaccard distance between a reference and itself must be zero.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual(0.0, dm.jaccardDistance("ref1", "ref1"))

    def testJaccardDistanceToIdentical(self):
        """
        The Jaccard distance between a reference and another with the same set
        of matching queries must be zero.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual(0.0, dm.jaccardDistance("ref1", "ref1"))

    def testJaccardDistanceWithNoQueriesInCommon(self):
        """
        The Jaccarddistance between two references that have no matching
        queries in common must be 1.0.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual(1.0, dm.jaccardDistance("ref1", "ref2"))

    def testJaccardDistanceWithOneQueryInCommon(self):
        """
        The Jaccard similarity between two references with one query in common
        is one over the number of queries that match them in total (four),
        i.e., 1/4 and the Jaccard distance is 1.0 minus this, or 3/4.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query4 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual(0.75, dm.jaccardDistance("ref1", "ref2"))

    def testJaccardDistanceWithTwoQueriesInCommon(self):
        """
        The Jaccard similarity between two references with two queries in
        common is two over the number of queries that match them in total
        (five), i.e., 2/5 and the Jaccard distance is 1.0 minus this, or 3/5.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query4 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query5 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query5 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual(0.6, dm.jaccardDistance("ref1", "ref2"))

    def testSoergelDistanceWithNegativeScore(self):
        """
        Soergel distance cannot be computed if a negative score is present.
        A ValueError must be raised in such cases.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:-50",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            error = (
                rf"^Alignment 1 in {filename!r} has tag 'AS' with "
                rf"negative value \(-50\)\.$"
            )
            self.assertRaisesRegex(
                ValueError, error, dm.addFile, filename, scoreTag="AS"
            )

    def testSoergelDistanceWithOneQueryInCommonNoScoreTag(self):
        """
        The Soergel similarity between two references with one query in common
        if no score tag was given is one over the number of queries that match
        them in total (four), i.e., 1/4 and the distance is 1.0 minus this, or
        3/4.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
                "query4 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename)

        self.assertEqual(0.75, dm.soergelDistance("ref1", "ref2"))

    def testSoergelDistanceWithNoQueryInCommon(self):
        """
        The Soergel similarity between two references with no queries in common
        when using a score tag given is the sum of the minimum scores (all are
        zero) over the sum of the maximum scores (50 + 10 + 60 + 30 = 150),
        i.e., zero, and the distance is 1.0 minus this, or 1.0.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:50",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
                "query4 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:30",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(1.0, dm.soergelDistance("ref1", "ref2"))

    def testSoergelDistanceToIdentical(self):
        """
        The Soergel similarity between two references with two queries in
        common with the same scores must be zero.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query2 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:20",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:20",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(0.0, dm.soergelDistance("ref1", "ref2"))

    def testSoergelDistanceSameQueriesDifferentScores(self):
        """
        The Soergel similarity between two references with two queries in
        common but with different scores is the sum of the minimum scores
        (10 + 15 = 25) over the sum of the maximum scores (30 + 70 = 100),
        or 1/4, and the distance is 1.0 minus this, or 3/4. The unrelated
        query3 and ref3 are ignored.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "@SQ SN:ref3 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:30",
                "query2 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:15",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:70",
                "query3 0 ref3 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:70",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(0.75, dm.soergelDistance("ref1", "ref2"))

    def testSoergelDistanceWithOneQueryInCommon(self):
        """
        The Soergel similarity between two references with one query in common
        when using a score tag given is the sum of the minimum scores (30) over
        the sum of the maximum scores (50 + 10 + 60 = 120), or 1/4, and the
        distance is 1.0 minus this, or 3/4.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:30",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:50",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(0.75, dm.soergelDistance("ref1", "ref2"))

    def testSoergelDistanceWithTwoQueriesInCommon(self):
        """
        The Soergel similarity between two references with two queries in
        common when using a score tag given is the sum of the minimum scores
        (10 + 20) over the sum of the maximum scores (50 + 10 + 60 = 120),
        or 1/4, and the distance is 1.0 minus this, or 3/4.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:50",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query3 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:20",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        self.assertEqual(0.75, dm.soergelDistance("ref1", "ref2"))

    def testSave(self):
        """
        The save method must write out the correct JSON.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:f:10.0",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:f:11.0",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:f:12.0",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        fp = StringIO()
        dm.save(fp)

        self.assertEqual(
            {
                "ref1": {
                    "query1": 10.0,
                },
                "ref2": {
                    "query1": 11.0,
                    "query2": 12.0,
                },
            },
            loads(fp.getvalue()),
        )

    def testLoad(self):
        """
        The load method must read the JSON and store it correctly.
        """
        data = {
            "ref1": {
                "query1": 10.0,
            },
            "ref2": {
                "query1": 11.0,
                "query2": 12.0,
            },
        }
        dm = DistanceMatrix()
        fp = StringIO(dumps(data))
        dm.load(fp)
        self.assertEqual(data, dm.scores)

    def testJaccardMatrixWithOneQueryInCommon(self):
        """
        The Jaccard similarity between two references is the number of reads in
        common (1) over the number of reads in the union (4) or 0.25. Check
        that the distance and similarity matrice have the right values.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:30",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:50",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
                "query4 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        # Test distance.
        matrix = dm.matrix(metric="jaccard")

        self.assertTrue(
            np.array_equal(
                [
                    [0.00, 0.75],
                    [0.75, 0.00],
                ],
                matrix,
            )
        )

        # Test similarity.
        matrix = dm.matrix(metric="jaccard", similarity=True)

        self.assertTrue(
            np.array_equal(
                [
                    [1.00, 0.25],
                    [0.25, 1.00],
                ],
                matrix,
            )
        )

    def testJaccardMatrixWithOneQueryInCommonReturnDict(self):
        """
        The Jaccard similarity between two references is the number of reads in
        common (1) over the number of reads in the union (4) or 0.25. Check
        that the distance and similarity matrice have the right values, asking
        for a dictionary to be returned.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:30",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:50",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
                "query4 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        # Test distance.
        matrix = dm.matrix(metric="jaccard", returnDict=True)

        self.assertEqual(
            {
                "ref1": {
                    "ref1": 0.00,
                    "ref2": 0.75,
                },
                "ref2": {
                    "ref1": 0.75,
                    "ref2": 0.00,
                },
            },
            matrix,
        )

        # Test similarity.
        matrix = dm.matrix(metric="jaccard", similarity=True, returnDict=True)

        self.assertEqual(
            {
                "ref1": {
                    "ref1": 1.00,
                    "ref2": 0.25,
                },
                "ref2": {
                    "ref1": 0.25,
                    "ref2": 1.00,
                },
            },
            matrix,
        )

    def testSoergelMatrixWithOneQueryInCommon(self):
        """
        The Soergel similarity between two references with one query in common
        when using a score tag given is the sum of the minimum scores (30) over
        the sum of the maximum scores (50 + 10 + 60 = 120), or 1/4, and the
        distance is 1.0 minus this, or 3/4. Check that the distance and
        similarity matrice have the right values.
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:30",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:50",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        # Test distance.
        matrix = dm.matrix()

        self.assertTrue(
            np.array_equal(
                [
                    [0.00, 0.75],
                    [0.75, 0.00],
                ],
                matrix,
            )
        )

        # Test similarity.
        matrix = dm.matrix(similarity=True)

        self.assertTrue(
            np.array_equal(
                [
                    [1.00, 0.25],
                    [0.25, 1.00],
                ],
                matrix,
            )
        )

    def testSoergelMatrixWithOneQueryInCommonExplicitReferenceIds(self):
        """
        The Soergel similarity between two references with one query in common
        when using a score tag given is the sum of the minimum scores (30) over
        the sum of the maximum scores (50 + 10 + 60 = 120), or 1/4, and the
        distance is 1.0 minus this, or 3/4. Check that the distance and
        similarity matrice have the right values when an explicit list of
        reference ids is passed (other references must be ignored).
        """
        data = "\n".join(
            [
                "@SQ SN:ref1 LN:10",
                "@SQ SN:ref2 LN:10",
                "@SQ SN:ref3 LN:10",
                "query1 0 ref1 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:30",
                "query1 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:50",
                "query2 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:10",
                "query3 0 ref2 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
                "query3 0 ref3 2 60 2=2X2M * 0 0 TCTAGG 123456 AS:i:60",
            ]
        ).replace(" ", "\t")

        dm = DistanceMatrix()
        with dataFile(data) as filename:
            dm.addFile(filename, scoreTag="AS")

        # Test distance.
        matrix = dm.matrix(referenceIds=("ref1", "ref2"))

        self.assertTrue(
            np.array_equal(
                [
                    [0.00, 0.75],
                    [0.75, 0.00],
                ],
                matrix,
            )
        )

        # Test similarity.
        matrix = dm.matrix(referenceIds=("ref1", "ref2"), similarity=True)

        self.assertTrue(
            np.array_equal(
                [
                    [1.00, 0.25],
                    [0.25, 1.00],
                ],
                matrix,
            )
        )


@dataclass
class AlignedRead:
    """
    Hold aligned read attributes that correspond to the pysam.AlignedRead class.

    We cannot just make an instance of that class using pure Python because
    it is implemented in C using data read from a BAM/SAM file and its attributes
    are read-only.
    """

    is_duplicate: bool = False
    is_qcfail: bool = False
    is_secondary: bool = False
    is_supplementary: bool = False
    query_name: str = ""
    reference_length: int = 0
    reference_start: int = 0


class TestReferenceReads(TestCase):
    """
    Test the ReferenceReads class.
    """

    def testAttributes(self):
        "A ReferenceReads instance must have attributes with the expected values."
        rr = ReferenceReads("id", 500)
        self.assertEqual("id", rr.id_)
        self.assertEqual(500, rr.length)
        # The read id sets must all be empty.
        for attr in (
            "coveredOffsets",
            "duplicate",
            "nonDuplicate",
            "primary",
            "qcFail",
            "readIds",
            "secondary",
            "supplementary",
        ):
            self.assertEqual(set(), getattr(rr, attr))

    def testQueryName(self):
        "Query names (read ids) must be added correctly."
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1"))
        rr.add(AlignedRead(query_name="q2"))
        self.assertEqual({"q1", "q2"}, rr.readIds)

    def testPrimary(self):
        "Primary (non-secondary, non-supplementary) reads must be added correctly."
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1"))
        self.assertEqual({"q1"}, rr.primary)

    def testSecondary(self):
        "Secondary reads must be added correctly."
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1", is_secondary=True))
        self.assertEqual({"q1"}, rr.secondary)

    def testSupplementar(self):
        "Supplementary reads must be added correctly."
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1", is_supplementary=True))
        self.assertEqual({"q1"}, rr.supplementary)

    def testDuplicate(self):
        "Duplicate reads must be added correctly."
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1", is_duplicate=True))
        self.assertEqual({"q1"}, rr.duplicate)

    def testNonDuplicate(self):
        "Non-duplicate reads must be added correctly."
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1"))
        self.assertEqual({"q1"}, rr.nonDuplicate)

    def testDuplicateId(self):
        "If a read id is repeated, the read id must be placed in the duplicate set."
        # See the comment about this possibility in the 'add' method of the
        # ReferenceReads class in ../dark/sam.py
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1"))
        rr.add(AlignedRead(query_name="q1"))
        self.assertEqual({"q1"}, rr.duplicate)
        self.assertEqual(set(), rr.nonDuplicate)

    def testQcFail(self):
        "QC failed reads must be added correctly."
        rr = ReferenceReads("id", 500)
        rr.add(AlignedRead(query_name="q1", is_qcfail=True))
        self.assertEqual({"q1"}, rr.qcFail)

    def testZeroCoverage(self):
        "If no reads are added, the coverage fraction must be 0.0"
        rr = ReferenceReads("id", 500)
        self.assertEqual(0.0, rr.coverage())

    def testCoverageTwentyPercent(self):
        "If one-fifth of the reference genome is covered, coverage() must return 0.2"
        rr = ReferenceReads("id", 500)
        read = AlignedRead(
            reference_start=100,
            reference_length=100,
        )
        rr.add(read)
        self.assertEqual(set(range(100, 200)), rr.coveredOffsets)
        self.assertEqual(0.2, rr.coverage())

    def testOverlappingCoverageEightyPercent(self):
        "If 80% of the reference genome is covered, coverage() must return 0.8"
        rr = ReferenceReads("id", 500)
        for count, (start, length) in enumerate(
            (
                (50, 100),
                (150, 100),
                (250, 100),
                (350, 100),
            )
        ):
            rr.add(
                AlignedRead(
                    query_name=str(count),
                    reference_start=start,
                    reference_length=length,
                )
            ),

        self.assertEqual(set(range(50, 450)), rr.coveredOffsets)
        self.assertEqual(0.8, rr.coverage())

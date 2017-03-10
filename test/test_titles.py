# Note: Tests for the TitlesAlignments class are in blast/test_titles.py
#       because that class needs a concrete (iterable)
#       dark.alignments.ReadsAlignments class passed to its __init__.  The
#       tests below test the simpler dark.titles classes, TitleAlignment
#       and TitleAlignments.

from collections import Counter
import six
import warnings
import platform
from unittest import TestCase

from dark.titles import TitleAlignment, TitleAlignments
from dark.reads import Read
from dark.hsp import HSP, LSP

_pypy = platform.python_implementation() == 'PyPy'


class WarningTestMixin(object):
    """
    Provide an assertion test which checks to see that a specified warning
    was raised.
    """
    # Taken from
    # http://stackoverflow.com/questions/3892218/
    # how-to-test-with-pythons-unittest-that-a-warning-has-been-thrown

    def assertWarns(self, warning, callable, *args, **kwds):
        with warnings.catch_warnings(record=True) as warning_list:
            warnings.simplefilter('always')
            callable(*args, **kwds)
            self.assertTrue(
                any(item.category == warning for item in warning_list))


class TestTitleAlignment(TestCase):
    """
    Test the TitleAlignment class.
    """

    def testExpectedAttributes(self):
        """
        An instance of TitleAlignment must have the expected attributes.
        """
        read = Read('id', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        self.assertEqual(read, titleAlignment.read)
        self.assertEqual([], titleAlignment.hsps)


class TestTitleAlignments(WarningTestMixin, TestCase):
    """
    Test the TitleAlignments class.
    """

    def testExpectedAttributes(self):
        """
        An instance of TitleAlignments must have the expected attributes.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        self.assertEqual('subject title', titleAlignments.subjectTitle)
        self.assertEqual(55, titleAlignments.subjectLength)
        self.assertEqual([], titleAlignments)

    def testAddAlignment(self):
        """
        It must be possible to add an alignment to an instance of
        TitleAlignments.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(read, titleAlignments[0].read)
        self.assertEqual([], titleAlignments[0].hsps)

    def testHSPs(self):
        """
        The hsps function must produce a list of all HSPs.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(14)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual([7, 14, 21],
                         [hsp.score.score for hsp in titleAlignments.hsps()])

    def testReadsEmpty(self):
        """
        The reads function must return an empty Reads instance if there are no
        reads for the title.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        self.assertEqual(0, len(list(titleAlignments.reads())))

    def testReads(self):
        """
        The reads function must return a Reads instance with the reads for
        the title.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(14)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read1 = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read1, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read2 = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read2, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual([read1, read2], list(titleAlignments.reads()))

    def testReadCountZero(self):
        """
        The readCount function must return zero if no reads matched a title.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        self.assertEqual(0, titleAlignments.readCount())

    def testReadCount(self):
        """
        The readCount function must indicate how many reads matched a title.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(14)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(2, titleAlignments.readCount())

    def testHspCountZero(self):
        """
        The hspCount function must return zero if no reads matched a title.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        self.assertEqual(0, titleAlignments.hspCount())

    def testHspCount(self):
        """
        The hspCount function must indicate how many HSPs were found in
        total for all the alignments to a title.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(14)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(3, titleAlignments.hspCount())

    def testMedianScoreWithNoAlignments(self):
        """
        The medianScore function must raise IndexError (due to no inputs)
        if there are no alignments matching a title.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        error = '^arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error,
                              titleAlignments.medianScore)

    def testMedianScoreWithNoHsps(self):
        """
        The medianScore function must raise ValueError if there are no HSPs.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        titleAlignments.addAlignment(titleAlignment)
        error = '^arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error,
                              titleAlignments.medianScore)

    def testMedianScoreOfTwo(self):
        """
        The medianScore function must return the median score for the HSPs in
        all the alignments matching a title when given 2 scores.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(15)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(11, titleAlignments.medianScore())

    def testMedianScoreOfThree(self):
        """
        The medianScore function must return the median score for the HSPs in
        all the alignments matching a title when given 3 scores.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(15)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(15, titleAlignments.medianScore())

    def testBestHspWithNoHsps(self):
        """
        The bestHsp function must raise ValueError if there are no HSPs.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        titleAlignments.addAlignment(titleAlignment)
        if _pypy:
            error = '^arg is an empty sequence$'
        else:
            error = '^max\(\) arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error, titleAlignments.bestHsp)

    def testBestHsp(self):
        """
        The bestHsp function must return the HSP with the best score for all
        the HSPs for all the alignments matching a title.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(15)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(hsp3, titleAlignments.bestHsp())

    def testWorstHspWithNoHsps(self):
        """
        The worstHsp function must raise ValueError if there are no HSPs.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        titleAlignments.addAlignment(titleAlignment)
        if _pypy:
            error = '^arg is an empty sequence$'
        else:
            error = '^min\(\) arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error,
                              titleAlignments.worstHsp)

    def testWorstHsp(self):
        """
        The worstHsp function must return the HSP with the worst score for all
        the HSPs for all the alignments matching a title.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(15)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(hsp1, titleAlignments.worstHsp())

    def testBetterThanFalse(self):
        """
        The hasScoreBetterThan function must return False if there is no HSP
        with a score better than the passed value.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(15)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertFalse(titleAlignments.hasScoreBetterThan(21))

    def testBetterThanTrue(self):
        """
        The hasScoreBetterThan function must return True if there is an HSP
        with a score better than the passed value.
        """
        hsp1 = HSP(7)
        hsp2 = HSP(15)
        hsp3 = HSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertTrue(titleAlignments.hasScoreBetterThan(20))

    def testCoverageNoReads(self):
        """
        The coverage method must return zero when a title alignments has no
        alignments (and therefore no coverage).
        """
        titleAlignments = TitleAlignments('subject title', 100)
        self.assertEqual(0.0, titleAlignments.coverage())

    def testFullCoverage(self):
        """
        The coverage method must return the correct value when the title is
        fully covered by its reads.
        """
        hsp1 = HSP(7, subjectStart=0, subjectEnd=50)
        hsp2 = HSP(8, subjectStart=50, subjectEnd=100)
        titleAlignments = TitleAlignments('subject title', 100)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(1.0, titleAlignments.coverage())

    def testPartialCoverage(self):
        """
        The coverage method must return the correct value when the title is
        partially covered by its reads.
        """
        hsp1 = HSP(7, subjectStart=10, subjectEnd=20)
        hsp2 = HSP(15, subjectStart=30, subjectEnd=40)
        hsp3 = HSP(21, subjectStart=50, subjectEnd=60)
        titleAlignments = TitleAlignments('subject title', 100)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(0.3, titleAlignments.coverage())

    def testFullCoverageCounts(self):
        """
        The coverageCounts method must return the correct result when the title
        is fully covered by its reads.
        """
        hsp1 = HSP(7, subjectStart=0, subjectEnd=5)
        hsp2 = HSP(8, subjectStart=5, subjectEnd=10)
        titleAlignments = TitleAlignments('subject title', 10)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        c = Counter([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(c, titleAlignments.coverageCounts())

    def testCoverageCounts(self):
        """
        The coverageCounts method must return the correct results when the
        title is partially covered by its reads.
        """
        hsp1 = HSP(7, subjectStart=1, subjectEnd=2)
        hsp2 = HSP(15, subjectStart=3, subjectEnd=4)
        hsp3 = HSP(21, subjectStart=5, subjectEnd=6)
        titleAlignments = TitleAlignments('subject title', 10)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        c = Counter([1, 3, 5])
        self.assertEqual(c, titleAlignments.coverageCounts())

    def testCoverageCountsOverlap(self):
        """
        The coverageCounts method must return the correct results when the
        title is partially covered by its reads that overlap.
        """
        hsp1 = HSP(7, subjectStart=1, subjectEnd=2)
        hsp2 = HSP(15, subjectStart=3, subjectEnd=6)
        hsp3 = HSP(21, subjectStart=5, subjectEnd=6)
        titleAlignments = TitleAlignments('subject title', 10)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        c = Counter([1, 3, 4, 5, 5])
        self.assertEqual(c, titleAlignments.coverageCounts())

    def testCoverageInfoNoReads(self):
        """
        When a title has no reads aligned to it, the coverageInfo method
        must return an empty result.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        coverage = titleAlignments.coverageInfo()
        self.assertEqual({}, coverage)

    def testCoverageInfoOneReadWithOneHSP(self):
        """
        When a title has one read with one HSP aligned to it, the coverageInfo
        method must return just the indices and bases from that read.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        hsp = HSP(15, subjectStart=3, subjectEnd=6, readMatchedSequence='CGT')
        read = Read('id1', 'AAACGT')
        titleAlignment = TitleAlignment(read, [hsp])
        titleAlignments.addAlignment(titleAlignment)
        coverage = titleAlignments.coverageInfo()
        self.assertEqual(
            {
                3: [(15, 'C')],
                4: [(15, 'G')],
                5: [(15, 'T')],
            },
            coverage)

    def testCoverageInfoOneReadWithTwoHSPs(self):
        """
        When a title has one read with two HSPs aligned to it, the coverageInfo
        method must return the correct indices and bases from that read.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        hsp1 = HSP(15, subjectStart=1, subjectEnd=4, readMatchedSequence='A-A')
        hsp2 = HSP(10, subjectStart=3, subjectEnd=6, readMatchedSequence='CGT')
        read = Read('id1', 'AAACGT')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        coverage = titleAlignments.coverageInfo()
        self.assertEqual(
            {
                1: [(15, 'A')],
                2: [(15, '-')],
                3: [(15, 'A'), (10, 'C')],
                4: [(10, 'G')],
                5: [(10, 'T')],
            },
            coverage)

    def testCoverageInfoTwoReadsWithThreeHSPs(self):
        """
        When a title has two reads (one with two HSPs, one with one) aligned
        to it, the coverageInfo method must return the correct indices and
        bases from the read.
        """
        titleAlignments = TitleAlignments('subject title', 55)

        # First read.
        hsp1 = HSP(15, subjectStart=1, subjectEnd=4, readMatchedSequence='A-A')
        hsp2 = HSP(10, subjectStart=3, subjectEnd=6, readMatchedSequence='CGT')
        read = Read('id1', 'AAACGT')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)

        # Second read.
        hsp1 = HSP(20, subjectStart=5, subjectEnd=10,
                   readMatchedSequence='CGGTA')
        read = Read('id2', 'AAACGTCGGTAAAA')
        titleAlignment = TitleAlignment(read, [hsp1])
        titleAlignments.addAlignment(titleAlignment)

        coverage = titleAlignments.coverageInfo()
        self.assertEqual(
            {
                1: [(15, 'A')],
                2: [(15, '-')],
                3: [(15, 'A'), (10, 'C')],
                4: [(10, 'G')],
                5: [(10, 'T'), (20, 'C')],
                6: [(20, 'G')],
                7: [(20, 'G')],
                8: [(20, 'T')],
                9: [(20, 'A')],
            },
            coverage)

    def testResidueCountsNoReads(self):
        """
        When a title has no reads aligned to it, the residueCounts method
        must return an empty result.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        counts = titleAlignments.residueCounts()
        self.assertEqual(0, len(counts))

    def testResidueCountsUnknownCaseConversion(self):
        """
        The residueCounts method must raise a ValueError when asked to do an
        unknown case conversion.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        error = "convertCaseTo must be one of 'none', 'lower', or 'upper'"
        six.assertRaisesRegex(
            self, ValueError, error, titleAlignments.residueCounts,
            convertCaseTo='xxx')

    def testResidueCountsOneReadOneHSP(self):
        """
        The residueCounts method must return the correct result when just one
        read with one HSP is aligned to a title.
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=0, readEnd=4, readStartInSubject=0,
                  readEndInSubject=4, subjectStart=0, subjectEnd=4,
                  readMatchedSequence='ACGT', subjectMatchedSequence='ACGT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                0: {'A': 1},
                1: {'C': 1},
                2: {'G': 1},
                3: {'T': 1},
            },
            titleAlignments.residueCounts())

    def testResidueCountsOneReadOneHSPPartialMatch(self):
        """
        The residueCounts method must return the correct result when just one
        read with one HSP is aligned to a title and only part of the read
        matched the subject (all the read bases are still counted and
        returned).
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=0, readEnd=2, readStartInSubject=0,
                  readEndInSubject=4, subjectStart=0, subjectEnd=4,
                  readMatchedSequence='ACGT', subjectMatchedSequence='ACGT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                0: {'A': 1},
                1: {'C': 1},
                2: {'G': 1},
                3: {'T': 1},
            },
            titleAlignments.residueCounts())

    def testResidueCountsOneReadTwoHSPsAtStartOfSubject(self):
        """
        The residueCounts method must return the correct result when just one
        read with two HSPs is aligned to a title and the leftmost HSP is
        aligned with the left edge of the subject.

        HSP1:       ACGT
        HSP2:        CGTT
        """
        read = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=0,
                   readEndInSubject=4, subjectStart=0, subjectEnd=4,
                   readMatchedSequence='ACGT', subjectMatchedSequence='ACGT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=1,
                   readEndInSubject=5, subjectStart=1, subjectEnd=5,
                   readMatchedSequence='CGTT', subjectMatchedSequence='CGTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                0: {'A': 1},
                1: {'C': 2},
                2: {'G': 2},
                3: {'T': 2},
                4: {'T': 1},
            },
            titleAlignments.residueCounts())

    def testResidueCountsOneReadTwoHSPsNotAtStartOfSubject(self):
        """
        The residueCounts method must return the correct result when just one
        read with two HSPs is aligned to a title and the leftmost HSP is not
        aligned with the left edge of the subject.

        HSP1:       ACGT
        HSP2:        CGTT
        """
        read = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=10,
                   readEndInSubject=14, subjectStart=10, subjectEnd=14,
                   readMatchedSequence='ACGT', subjectMatchedSequence='ACGT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=11,
                   readEndInSubject=15, subjectStart=11, subjectEnd=15,
                   readMatchedSequence='CGTT', subjectMatchedSequence='CGTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                10: {'A': 1},
                11: {'C': 2},
                12: {'G': 2},
                13: {'T': 2},
                14: {'T': 1},
            },
            titleAlignments.residueCounts())

    def testResidueCountsNoCaseConversion(self):
        """
        The residueCounts method must return the correct result when asked not
        to convert case.

        HSP1:       AcgT
        HSP2:        CGTT
        """
        read = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=10,
                   readEndInSubject=14, subjectStart=10, subjectEnd=14,
                   readMatchedSequence='AcgT', subjectMatchedSequence='ACGT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=11,
                   readEndInSubject=15, subjectStart=11, subjectEnd=15,
                   readMatchedSequence='CGTT', subjectMatchedSequence='CGTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                10: {'A': 1},
                11: {'C': 1, 'c': 1},
                12: {'G': 1, 'g': 1},
                13: {'T': 2},
                14: {'T': 1},
            },
            titleAlignments.residueCounts(convertCaseTo='none'))

    def testResidueCountsCaseConvertLower(self):
        """
        The residueCounts method must return the correct result when asked to
        convert residues to lower case.

        HSP1:       AcgT
        HSP2:        CGTT
        """
        read = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=10,
                   readEndInSubject=14, subjectStart=10, subjectEnd=14,
                   readMatchedSequence='AcgT', subjectMatchedSequence='ACGT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=11,
                   readEndInSubject=15, subjectStart=11, subjectEnd=15,
                   readMatchedSequence='CGTT', subjectMatchedSequence='CGTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                10: {'a': 1},
                11: {'c': 2},
                12: {'g': 2},
                13: {'t': 2},
                14: {'t': 1},
            },
            titleAlignments.residueCounts(convertCaseTo='lower'))

    def testResidueCountsCaseConvertUpper(self):
        """
        The residueCounts method must return the correct result when asked to
        convert residues to upper case.

        HSP1:       AcgT
        HSP2:        CGTT
        """
        read = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=10,
                   readEndInSubject=14, subjectStart=10, subjectEnd=14,
                   readMatchedSequence='AcgT', subjectMatchedSequence='ACGT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=11,
                   readEndInSubject=15, subjectStart=11, subjectEnd=15,
                   readMatchedSequence='CGTT', subjectMatchedSequence='CGTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                10: {'A': 1},
                11: {'C': 2},
                12: {'G': 2},
                13: {'T': 2},
                14: {'T': 1},
            },
            titleAlignments.residueCounts(convertCaseTo='upper'))

    def testResidueCountsCaseConvertUpperIsDefault(self):
        """
        The residueCounts method must convert to uppercase by default.

        HSP1:       AcgT
        HSP2:        CGTT
        """
        read = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=10,
                   readEndInSubject=14, subjectStart=10, subjectEnd=14,
                   readMatchedSequence='AcgT', subjectMatchedSequence='ACGT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=11,
                   readEndInSubject=15, subjectStart=11, subjectEnd=15,
                   readMatchedSequence='CGTT', subjectMatchedSequence='CGTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                10: {'A': 1},
                11: {'C': 2},
                12: {'G': 2},
                13: {'T': 2},
                14: {'T': 1},
            },
            titleAlignments.residueCounts())

    def testResidueCountsTwoReadsTwoHSPsLeftOverhang(self):
        """
        The residueCounts method must return the correct result when two
        reads, each with one HSP are aligned to a title and the leftmost HSP
        is aligned before the left edge of the subject (i.e, will include
        negative subject offsets).

        Subject:      GTT
        HSP1:       ACGT
        HSP2:        CGTT
        """
        read1 = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=-2,
                   readEndInSubject=2, subjectStart=0, subjectEnd=2,
                   readMatchedSequence='GT', subjectMatchedSequence='GT')
        read2 = Read('id', 'CGTT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=-1,
                   readEndInSubject=3, subjectStart=0, subjectEnd=3,
                   readMatchedSequence='GTT', subjectMatchedSequence='GTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read1, [hsp1])
        titleAlignments.addAlignment(titleAlignment)
        titleAlignment = TitleAlignment(read2, [hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                -2: {'A': 1},
                -1: {'C': 2},
                0: {'G': 2},
                1: {'T': 2},
                2: {'T': 1},
            },
            titleAlignments.residueCounts())

    def testResidueCountsOneReadTwoHSPsNotOverlapping(self):
        """
        The residueCounts method must return the correct result when just one
        read with two HSPs is aligned to a title and the HSPs do not overlap
        one another.

        HSP1:    ACGT
        HSP2:              CGTT
        """
        read = Read('id', 'ACGT')
        hsp1 = HSP(33, readStart=0, readEnd=4, readStartInSubject=0,
                   readEndInSubject=4, subjectStart=0, subjectEnd=4,
                   readMatchedSequence='ACGT', subjectMatchedSequence='ACGT')
        hsp2 = HSP(33, readStart=0, readEnd=4, readStartInSubject=10,
                   readEndInSubject=14, subjectStart=10, subjectEnd=14,
                   readMatchedSequence='CGTT', subjectMatchedSequence='CGTT')
        titleAlignments = TitleAlignments('subject title', 55)
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(
            {
                0: {'A': 1},
                1: {'C': 1},
                2: {'G': 1},
                3: {'T': 1},
                10: {'C': 1},
                11: {'G': 1},
                12: {'T': 1},
                13: {'T': 1},
            },
            titleAlignments.residueCounts())

    def testSummaryWhenEmpty(self):
        """
        If summary is called on an instance of TitleAlignments with no
        alignments a ValueError must be raised.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        error = '^max\(\) arg is an empty sequence$'
        six.assertRaisesRegex(self, ValueError, error, titleAlignments.summary)

    def testSummary(self):
        """
        The summary method must return the correct result.
        """
        titleAlignments = TitleAlignments('subject title', 10)
        titleAlignments.addAlignment(
            TitleAlignment(Read('id1', 'ACGT'), [
                HSP(30, subjectStart=0, subjectEnd=2),
            ]))
        titleAlignments.addAlignment(
            TitleAlignment(Read('id2', 'ACGT'), [
                HSP(55, subjectStart=2, subjectEnd=4),
                HSP(40, subjectStart=8, subjectEnd=9),
            ]))
        self.assertEqual(
            {
                'bestScore': 55,
                'coverage': 0.5,
                'hspCount': 3,
                'medianScore': 40,
                'readCount': 2,
                'subjectLength': 10,
                'subjectTitle': 'subject title',
            },
            titleAlignments.summary())


class TestTitleAlignmentsLSP(TestCase):
    """
    Test the TitleAlignments class using LSPs. The only tests here are ones
    that depend on lower scores being better.
    """

    def testBestHsp(self):
        """
        The bestHsp function must return the HSP with the best score for the
        HSPs all the alignments matching a title.
        """
        hsp1 = LSP(7)
        hsp2 = LSP(15)
        hsp3 = LSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(hsp1, titleAlignments.bestHsp())

    def testWorstHsp(self):
        """
        The worstHsp function must return the HSP with the worst score for all
        the HSPs for all the alignments matching a title.
        """
        hsp1 = LSP(7)
        hsp2 = LSP(15)
        hsp3 = LSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(hsp3, titleAlignments.worstHsp())

    def testBetterThanFalse(self):
        """
        The hasScoreBetterThan function must return False if there is no HSP
        with a score better than the passed value.
        """
        hsp1 = LSP(7)
        hsp2 = LSP(15)
        hsp3 = LSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertFalse(titleAlignments.hasScoreBetterThan(7))

    def testBetterThanTrue(self):
        """
        The hasScoreBetterThan function must return True if there is an HSP
        with a score better than the passed value.
        """
        hsp1 = LSP(7)
        hsp2 = LSP(15)
        hsp3 = LSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertTrue(titleAlignments.hasScoreBetterThan(9))

    def testReadIdsEmpty(self):
        """
        The readIds function must return the empty set if no reads matched a
        title.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        self.assertEqual(0, len(titleAlignments.readIds()))

    def testReadIds(self):
        """
        The readIds function must return the set of read ids for the alignments
        matching a title.
        """
        hsp1 = LSP(7)
        hsp2 = LSP(15)
        hsp3 = LSP(21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(set(['id1', 'id2']), titleAlignments.readIds())

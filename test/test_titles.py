# Note: Tests for the TitlesAlignments class are in blast/test_titles.py
#       because that class needs a concrete (iterable)
#       dark.alignments.ReadsAlignments class passed to its __init__.  The
#       tests below test the simpler dark.titles classes, TitleAlignment
#       and TitleAlignments.

import warnings
from unittest import TestCase

from dark.titles import TitleAlignment, TitleAlignments
from dark.reads import Read
from dark.hsp import HSP, LSP


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
        self.assertEqual(0, len(titleAlignments.reads()))

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
        The medianScore function must issue a warning (due to no inputs)
        if there are no alignments matching a title.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        self.assertWarns(RuntimeWarning, titleAlignments.medianScore)

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

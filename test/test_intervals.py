from unittest import TestCase

from dark.intervals import OffsetAdjuster, ReadIntervals


class TestReadIntervals(TestCase):
    """
    Test the ReadIntervals class
    """

    EMPTY = ReadIntervals.EMPTY
    FULL = ReadIntervals.FULL

    def testEmpty(self):
        """
        When no intervals are added, walk should just return one empty
        interval that spans the entire rangoe from 0 to the sequence
        length.
        """
        ri = ReadIntervals(100)
        self.assertEqual(
            [
                (self.EMPTY, (0, 100))
            ],
            list(ri.walk()))

    def testOneIntervalExactCovering(self):
        """
        If there is a single interval that spans the whole hit exactly, just
        that one interval should be returned by walk, and it should be full.
        """
        ri = ReadIntervals(100)
        ri.add(0, 100)
        self.assertEqual(
            [
                (self.FULL, (0, 100))
            ],
            list(ri.walk()))

    def testOneIntervalCoveringAllExtendingLeft(self):
        """
        If there is a single interval that spans the whole hit, including
        going negative to the left, that one interval should be returned by
        walk, and it should be full.
        """
        ri = ReadIntervals(100)
        ri.add(-10, 100)
        self.assertEqual(
            [
                (self.FULL, (-10, 100))
            ],
            list(ri.walk()))

    def testOneIntervalCoveringAllExtendingRight(self):
        """
        If there is a single interval that spans the whole hit, including
        going beyond the hit to the right, that one interval should be returned
        by walk, and it should be full.
        """
        ri = ReadIntervals(100)
        ri.add(0, 110)
        self.assertEqual(
            [
                (self.FULL, (0, 110))
            ],
            list(ri.walk()))

    def testOneIntervalCoveringAllExtendingBoth(self):
        """
        If there is a single interval that spans the whole hit, including
        starting before zero and also going beyond the hit to the right, that
        one interval should be returned by walk, and it should be full.
        """
        ri = ReadIntervals(100)
        ri.add(-10, 110)
        self.assertEqual(
            [
                (self.FULL, (-10, 110))
            ],
            list(ri.walk()))

    def testOneIntervalStartingAtZero(self):
        """
        If there is a single interval that starts at zero but doesn't
        cover the whole hit, we should get 2 intervals back from walk,
        a full one and then an empty.
        """
        ri = ReadIntervals(100)
        ri.add(0, 50)
        self.assertEqual(
            [
                (self.FULL, (0, 50)),
                (self.EMPTY, (50, 100)),
            ],
            list(ri.walk()))

    def testOneIntervalStartingBeforeZero(self):
        """
        If there is a single interval that starts before zero but doesn't
        cover the whole hit, we should get 2 intervals back from walk,
        a full one and then an empty.
        """
        ri = ReadIntervals(100)
        ri.add(-50, 50)
        self.assertEqual(
            [
                (self.FULL, (-50, 50)),
                (self.EMPTY, (50, 100)),
            ],
            list(ri.walk()))

    def testOneIntervalEndingAtHitEnd(self):
        """
        If there is a single interval that ends at the end of the hit
        but doesn't start at zero, we should get 2 intervals back from walk,
        an empty then a full.
        """
        ri = ReadIntervals(100)
        ri.add(50, 100)
        self.assertEqual(
            [
                (self.EMPTY, (0, 50)),
                (self.FULL, (50, 100)),
            ],
            list(ri.walk()))

    def testOneIntervalEndingAfterHitEnd(self):
        """
        If there is a single interval that ends after the end of the hit
        but doesn't start at zero, we should get 2 intervals back from walk,
        an empty then a full.
        """
        ri = ReadIntervals(100)
        ri.add(50, 150)
        self.assertEqual(
            [
                (self.EMPTY, (0, 50)),
                (self.FULL, (50, 150)),
            ],
            list(ri.walk()))

    def testOneIntervalInMiddle(self):
        """
        If there is a single interval in the middle of the hit, we
        should get 3 intervals back from walk: empty, full, empty.

        """
        ri = ReadIntervals(100)
        ri.add(50, 60)
        self.assertEqual(
            [
                (self.EMPTY, (0, 50)),
                (self.FULL, (50, 60)),
                (self.EMPTY, (60, 100)),
            ],
            list(ri.walk()))

    def testTwoOverlappingIntervalsInMiddle(self):
        """
        If there are two overlapping intervals in the middle of the hit, we
        should get 3 intervals back from walk: empty, full, empty.
        """
        ri = ReadIntervals(100)
        ri.add(50, 60)
        ri.add(55, 70)
        self.assertEqual(
            [
                (self.EMPTY, (0, 50)),
                (self.FULL, (50, 70)),
                (self.EMPTY, (70, 100)),
            ],
            list(ri.walk()))

    def testThreeOverlappingIntervalsInMiddle(self):
        """
        If there are three overlapping intervals in the middle of the hit, we
        should get 3 intervals back from walk: empty, full, empty.
        """
        ri = ReadIntervals(100)
        ri.add(50, 60)
        ri.add(65, 75)
        ri.add(55, 70)
        self.assertEqual(
            [
                (self.EMPTY, (0, 50)),
                (self.FULL, (50, 75)),
                (self.EMPTY, (75, 100)),
            ],
            list(ri.walk()))

    def testPairOfTwoOverlappingIntervals(self):
        """
        If there are two sets of two overlapping intervals in the middle of
        the hit, we should get 5 intervals back from walk.
        """
        ri = ReadIntervals(100)
        # First overlapping pair, 50-70.
        ri.add(50, 60)
        ri.add(55, 70)
        # First overlapping pair, 80-95.
        ri.add(80, 90)
        ri.add(85, 95)
        self.assertEqual(
            [
                (self.EMPTY, (0, 50)),
                (self.FULL, (50, 70)),
                (self.EMPTY, (70, 80)),
                (self.FULL, (80, 95)),
                (self.EMPTY, (95, 100)),
            ],
            list(ri.walk()))

    def testOverlappingIntervalsThatCoverEverything(self):
        """
        If there are sets of overlapping intervals that cover the whole hit,
        we should get 1 full interval back from walk.
        """
        ri = ReadIntervals(100)
        ri.add(-10, 20)
        ri.add(15, 40)
        ri.add(40, 70)
        ri.add(66, 89)
        ri.add(77, 93)
        ri.add(70, 110)
        self.assertEqual(
            [
                (self.FULL, (-10, 110))
            ],
            list(ri.walk()))


class TestOffsetAdjuster(TestCase):
    """
    Tests for the OffsetAdjuster class.
    """

    def testEmpty(self):
        """
        When no intervals are given, the reduction should be for the full hit
        width.
        """
        adjuster = OffsetAdjuster()
        self.assertEqual([], adjuster.adjustments())
        self.assertEqual(0, adjuster.adjustOffset(0))

    def testNoReads(self):
        """
        When no reads are added to an interval, the reduction should be for
        the full hit width.
        """
        ri = ReadIntervals(64)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
                (64, 58)
            ],
            adjuster.adjustments())
        self.assertEqual(6, adjuster.adjustOffset(64))

    def testOneReadThatExactlyCoversHit(self):
        """
        When one read is given that exactly covers the hit, there should
        be no length reductions.
        """
        ri = ReadIntervals(106)
        ri.add(0, 106)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
            ],
            adjuster.adjustments())
        self.assertEqual(106, adjuster.adjustOffset(106))

    def testOneReadThatExceedsHitOnBothEnds(self):
        """
        When one read is given that exceeds the hit at both ends, there should
        be no length reductions.
        """
        ri = ReadIntervals(106)
        ri.add(-100, 200)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
            ],
            adjuster.adjustments())
        self.assertEqual(106, adjuster.adjustOffset(106))

    def testOneReadAtStart(self):
        """
        When one read is added to the start of an interval, there should be one
        reduction for the empty section after the read.
        """
        ri = ReadIntervals(228)
        ri.add(0, 100)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
                (228, 121),
            ],
            adjuster.adjustments())
        self.assertEqual(107, adjuster.adjustOffset(228))

    def testOneReadBeforeStart(self):
        """
        When one read is added to the start of an interval before zero, there
        should be one reduction for the empty section after the read.
        """
        ri = ReadIntervals(228)
        ri.add(-10, 100)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
                (228, 121),
            ],
            adjuster.adjustments())
        self.assertEqual(107, adjuster.adjustOffset(228))

    def testOneReadAtEnd(self):
        """
        When one read is added to the end of an interval, there should be one
        reduction for the empty section before the read.
        """
        ri = ReadIntervals(228)
        ri.add(128, 228)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
                (128, 121),
            ],
            adjuster.adjustments())
        self.assertEqual(107, adjuster.adjustOffset(228))

    def testOneReadAfterEnd(self):
        """
        When one read is added to the end of an interval, going beyond the end
        of the hit, there should be one reduction for the empty section before
        the read.
        """
        ri = ReadIntervals(228)
        ri.add(128, 250)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
                (128, 121),
            ],
            adjuster.adjustments())
        self.assertEqual(107, adjuster.adjustOffset(228))

    def testOneReadInMiddle(self):
        """
        When one read is added to the middle of an interval, there should be
        two reductions.
        """
        ri = ReadIntervals(106)
        ri.add(32, 42)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
                (32, 27),
                (106, 58),
            ],
            adjuster.adjustments())
        self.assertEqual(106 - 27 - 58, adjuster.adjustOffset(106))

    def testTwoReadsInMiddle(self):
        """
        When two reads are added to the middle of an interval, there should be
        three reductions (after first empty area, after 2nd empty area, after
        final empty area.
        """
        ri = ReadIntervals(132)
        ri.add(32, 42)
        ri.add(58, 68)
        adjuster = OffsetAdjuster(ri)
        self.assertEqual(
            [
                (32, 27),
                (58, 12),
                (132, 58),
            ],
            adjuster.adjustments())
        self.assertEqual(132 - 27 - 12 - 58, adjuster.adjustOffset(132))

        # Test an HSP at the beginning is unchanged.
        self.assertEqual(
            {
                'queryEnd': 10,
                'queryStart': 0,
                'subjectEnd': 10,
                'subjectStart': 0,
            },
            adjuster.adjustNormalizedHSP({
                'queryEnd': 10,
                'queryStart': 0,
                'subjectEnd': 10,
                'subjectStart': 0,
            }))

        # Test an HSP in the first read region.
        self.assertEqual(
            {
                'queryEnd': 15,
                'queryStart': 5,
                'subjectEnd': 13,
                'subjectStart': 8,
            },
            adjuster.adjustNormalizedHSP({
                'queryEnd': 42,
                'queryStart': 32,
                'subjectEnd': 40,
                'subjectStart': 35,
            }))

        # Test an HSP in the second read region.
        self.assertEqual(
            {
                'queryEnd': 29,
                'queryStart': 19,
                'subjectEnd': 27,
                'subjectStart': 21,
            },
            adjuster.adjustNormalizedHSP({
                'queryEnd': 68,
                'queryStart': 58,
                'subjectEnd': 66,
                'subjectStart': 60,
            }))

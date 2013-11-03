from unittest import TestCase

from dark.intervals import ReadIntervals


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

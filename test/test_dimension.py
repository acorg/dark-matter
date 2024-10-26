from unittest import TestCase

from dark.dimension import dimensionalIterator


class TestDimensionalIterator(TestCase):
    """
    Test the dimensionalIterator generator.
    """

    def testNoDimensions(self):
        self.assertEqual((), tuple(dimensionalIterator(())))

    def testZeroDimension(self):
        """
        Passing a zero dimension must result in a ValueError with the expected
        error string.
        """
        error = r"^Dimensions not all positive! \(2, 0, 3\)$"
        self.assertRaisesRegex(ValueError, error, next, dimensionalIterator((2, 0, 3)))

    def testNegativeDimension(self):
        """
        Passing a negative dimension must result in a ValueError with the
        expected error string.
        """
        error = r"^Dimensions not all positive! \(2, -1, 3\)$"
        self.assertRaisesRegex(ValueError, error, next, dimensionalIterator((2, -1, 3)))

    def testLimitedTo0Items(self):
        self.assertEqual((), tuple(dimensionalIterator((2, 2), maxItems=0)))

    def test1X2(self):
        self.assertEqual(((0, 0), (0, 1)), tuple(dimensionalIterator((1, 2))))

    def test2X2(self):
        self.assertEqual(
            ((0, 0), (0, 1), (1, 0), (1, 1)), tuple(dimensionalIterator((2, 2)))
        )

    def test2X3WithStartValueTooShort(self):
        """
        If passed a too-short start value, a ValueError must be raised.
        """
        error = (
            r"^The start list is of length 2, but the number of dimensions "
            r"\(3\) differs\.$"
        )
        self.assertRaisesRegex(
            ValueError, error, next, dimensionalIterator((2, 3, 4), start=(0, 0))
        )

    def test2X3WithStartValueTooLong(self):
        """
        If passed a too-long start value, a ValueError must be raised.
        """
        error = (
            r"^The start list is of length 4, but the number of dimensions "
            r"\(3\) differs\.$"
        )
        self.assertRaisesRegex(
            ValueError, error, next, dimensionalIterator((2, 3, 4), start=(0, 0, 0, 0))
        )

    def test2X3WithStartValue(self):
        self.assertEqual(
            ((0, 2), (1, 0), (1, 1), (1, 2)),
            tuple(dimensionalIterator((2, 3), start=(0, 2))),
        )

    def test2X3WithStartValueAndLimit(self):
        self.assertEqual(
            ((0, 2), (1, 0), (1, 1)),
            tuple(dimensionalIterator((2, 3), start=(0, 2), maxItems=3)),
        )

    def test2X2LimitedTo4Items(self):
        """
        If the passed maximum number of items is the same as the number
        of items we'd expect to get with no limit, we should get all
        items.
        """
        self.assertEqual(
            ((0, 0), (0, 1), (1, 0), (1, 1)),
            tuple(dimensionalIterator((2, 2), maxItems=4)),
        )

    def test2X2LimitedTo3Items(self):
        self.assertEqual(
            ((0, 0), (0, 1), (1, 0)), tuple(dimensionalIterator((2, 2), maxItems=3))
        )

    def testStarX2LimitedTo10Items(self):
        self.assertEqual(
            (
                (0, 0),
                (0, 1),
                (1, 0),
                (1, 1),
                (2, 0),
                (2, 1),
                (3, 0),
                (3, 1),
                (4, 0),
                (4, 1),
            ),
            tuple(dimensionalIterator(("*", 2), maxItems=10)),
        )

    def testStarX3LimitedTo10Items(self):
        self.assertEqual(
            (
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
                (1, 2),
                (2, 0),
                (2, 1),
                (2, 2),
                (3, 0),
            ),
            tuple(dimensionalIterator(("*", 3), maxItems=10)),
        )

    def test3X2XStarX4(self):
        """
        If '*' is used in a dimension that's not the first, iteration
        should get 'stuck' there. I.e., the earlier dimensions will always
        return zero as the '*' dimension is never exhauted.
        """
        self.assertEqual(
            (
                (0, 0, 0, 0),
                (0, 0, 0, 1),
                (0, 0, 0, 2),
                (0, 0, 0, 3),
                (0, 0, 1, 0),
                (0, 0, 1, 1),
                (0, 0, 1, 2),
                (0, 0, 1, 3),
                (0, 0, 2, 0),
                (0, 0, 2, 1),
                (0, 0, 2, 2),
                (0, 0, 2, 3),
                (0, 0, 3, 0),
                (0, 0, 3, 1),
                (0, 0, 3, 2),
                (0, 0, 3, 3),
                (0, 0, 4, 0),
                (0, 0, 4, 1),
                (0, 0, 4, 2),
                (0, 0, 4, 3),
            ),
            tuple(dimensionalIterator((3, 2, "*", 4), maxItems=20)),
        )

    def test2X3X4(self):
        self.assertEqual(
            (
                (0, 0, 0),
                (0, 0, 1),
                (0, 0, 2),
                (0, 0, 3),
                (0, 1, 0),
                (0, 1, 1),
                (0, 1, 2),
                (0, 1, 3),
                (0, 2, 0),
                (0, 2, 1),
                (0, 2, 2),
                (0, 2, 3),
                (1, 0, 0),
                (1, 0, 1),
                (1, 0, 2),
                (1, 0, 3),
                (1, 1, 0),
                (1, 1, 1),
                (1, 1, 2),
                (1, 1, 3),
                (1, 2, 0),
                (1, 2, 1),
                (1, 2, 2),
                (1, 2, 3),
            ),
            tuple(dimensionalIterator((2, 3, 4))),
        )

    def test1X1X1(self):
        self.assertEqual(((0, 0, 0),), tuple(dimensionalIterator((1, 1, 1))))

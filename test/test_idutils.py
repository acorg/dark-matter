from unittest import TestCase

from dark.idutils import ids


class TestIds(TestCase):
    """
    Test the ids function.
    """

    def testStartAndStopStringsOfUnequalLength(self):
        """
        If the passed start and end string are not of the same length, a ValueError
        must be raised.
        """
        error = r"^Range specifiers must be of equal length\.$"
        self.assertRaisesRegex(ValueError, error, tuple, ids("xxx", "yyyy"))

    def testNothing(self):
        """
        If no ids are requested, an empty list must result.
        """
        self.assertEqual((), tuple(ids("Pool-AB", maxResults=0)))

    def testMaxCount(self):
        """
        A non-zero maximum count must work as expected with just a start
        value passed.
        """
        self.assertEqual(("aa", "ab", "ac", "ad"), tuple(ids("aa", maxResults=4)))

    def testMaxCountCutsUnfinishedIds(self):
        """
        A non-zero maximum count must work as expected when the number of ids
        that would be returned exceeds the count.
        """
        self.assertEqual(("6", "7", "8"), tuple(ids(6, 9, maxResults=3)))

    def testTwoNakedIntsAscending(self):
        """
        Test two naked integers in ascending order.
        """
        self.assertEqual(("5", "6", "7"), tuple(ids(5, 7)))

    def testTwoNakedIntsDescending(self):
        """
        Test two naked integers in descending order.
        """
        self.assertEqual(("5", "6", "7"), tuple(ids(7, 5)))

    def testTwoIntsInAStringAscending(self):
        """
        Test two integers in ascending order in a string.
        """
        self.assertEqual(("id-5", "id-6", "id-7"), tuple(ids("id-5", "id-7")))

    def testTwoIntsInAStringDescending(self):
        """
        Test two integers given in descending order in a string.
        """
        self.assertEqual(("id-5", "id-6", "id-7"), tuple(ids("id-7", "id-5")))

    def testTwoIntsAscendingExplicitPrefixPassingStrings(self):
        """
        Test two integers in ascending order when arguments are passed as strings
        and an explicit prefix is given.
        """
        self.assertEqual(("id-5", "id-6", "id-7"), tuple(ids("5", "7", prefix="id-")))

    def testTwoIntsAscendingExplicitPrefixPassingInts(self):
        """
        Test two integers in ascending order when arguments are passed as ints
        and an explicit prefix is given.
        """
        self.assertEqual(("id-5", "id-6", "id-7"), tuple(ids(5, 7, prefix="id-")))

    def testTwoIntsDescendingExplicitPrefix(self):
        """
        Test two integers given in descending order when an explicit prefix is given.
        """
        self.assertEqual(("id-5", "id-6", "id-7"), tuple(ids(7, 5, prefix="id-")))

    def testIntermediateSymbols(self):
        """
        It must be possible to put intermediate (constant) symbols into ids.
        """
        self.assertEqual(
            ("id-5-7", "id-5-8", "id-5-9", "id-6-0", "id-6-1"),
            tuple(ids("id-5-7", "id-6-1")),
        )

    def testIntsWithNoLimit(self):
        """
        A starting int with no limit must iterate to the next power of ten.
        """
        self.assertEqual(("93", "94", "95", "96", "97", "98", "99"), tuple(ids(93)))

    def testTwoLetterPoolIdsWithStartAndNoStop(self):
        """
        Check that example fictitional pool ids work as expected when a start is
        given but a stop is not.
        """
        self.assertEqual(
            ("Pool-ZW", "Pool-ZX", "Pool-ZY", "Pool-ZZ"),
            tuple(ids("ZW", prefix="Pool-")),
        )

    def testTwoLetterPoolIdsWithStartAndStop(self):
        """
        Check that example fictitional pool ids work as expected when a start and
        stop id are both given.
        """
        self.assertEqual(
            ("Pool-BW", "Pool-BX", "Pool-BY", "Pool-BZ", "Pool-CA", "Pool-CB"),
            tuple(ids("Pool-BW", "Pool-CB")),
        )

    def testRandomThreeLetterStringsWithNoPrefix(self):
        """
        Check that an arbitrary series of three-letter strings can be generated.
        """
        self.assertEqual(
            (
                "hnx",
                "hny",
                "hnz",
                "hoa",
                "hob",
                "hoc",
                "hod",
                "hoe",
                "hof",
                "hog",
                "hoh",
                "hoi",
                "hoj",
                "hok",
                "hol",
                "hom",
                "hon",
                "hoo",
                "hop",
            ),
            tuple(ids("hnx", "hop", prefix="")),
        )

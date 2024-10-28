import bz2
import gzip
import builtins
from unittest import TestCase
from unittest.mock import mock_open, patch
from collections import Counter
from io import BytesIO, StringIO

from dark.utils import (
    asHandle,
    baseCountsToStr,
    countPrint,
    intsToIntervals,
    intsToStringIntervals,
    matchOffset,
    median,
    nucleotidesToStr,
    numericallySortFilenames,
    parseRangeExpression,
    parseRangeString,
    pct,
    take,
)


class TestNumericallySortFilenames(TestCase):
    """
    Test the numericallySortFilenames function.
    """

    def testNoNames(self):
        """
        An empty list must be returned when an empty list is given.
        """
        self.assertEqual([], numericallySortFilenames([]))

    def testOneNonNumericName(self):
        """
        A list with a single non-numeric name should result in that same
        name being returned.
        """
        self.assertEqual(["hey"], numericallySortFilenames(["hey"]))

    def testOneNumericName(self):
        """
        A list with a single numeric name should result in that same
        name being returned.
        """
        self.assertEqual(["3.json"], numericallySortFilenames(["3.json"]))

    def testSeveralNames(self):
        """
        A list with several numeric names should result in a correctly
        sorted list of names being returned.
        """
        self.assertEqual(
            ["1.json", "2.json", "3.json"],
            numericallySortFilenames(["3.json", "1.json", "2.json"]),
        )

    def testSeveralNamesWithUnequalPrefixLengths(self):
        """
        A list with several numeric names whose numeric prefixes differ
        in length should result in a correctly sorted list of names being
        returned.
        """
        self.assertEqual(
            ["2.json", "3.json", "21.json", "35.json", "250.json"],
            numericallySortFilenames(
                ["3.json", "21.json", "35.json", "250.json", "2.json"]
            ),
        )

    def testBasename(self):
        """
        Sorting must be according to file basename.
        """
        self.assertEqual(
            [
                "../output/2.json",
                "../output/3.json",
                "../output/21.json",
                "../output/35.json",
                "../output/250.json",
            ],
            numericallySortFilenames(
                [
                    "../output/3.json",
                    "../output/21.json",
                    "../output/35.json",
                    "../output/250.json",
                    "../output/2.json",
                ]
            ),
        )


class TestMedian(TestCase):
    """
    Tests for the median function.
    """

    def testEmptyArgRaises(self):
        """
        An empty list must cause median to raise ValueError.
        """
        error = "^arg is an empty sequence$"
        self.assertRaisesRegex(ValueError, error, median, [])

    def testMedianOfOne(self):
        """
        The median function must work on a list of length one.
        """
        self.assertEqual(3, median([3]))

    def testMedianOfTwo(self):
        """
        The median function must work on a list of length two.
        """
        self.assertEqual(4.5, median([3.1, 5.9]))

    def testMedianOfThree(self):
        """
        The median function must work on a list of length threee.
        """
        self.assertEqual(5.9, median([3.1, 7.6, 5.9]))

    def testMedianOfFour(self):
        """
        The median function must work on a list of length four.
        """
        self.assertEqual(4.5, median([3.1, 1.3, 7.6, 5.9]))

    def testMedianOfFive(self):
        """
        The median function must work on a list of length five.
        """
        self.assertEqual(5.9, median([3.1, 1.3, 7.6, 9.9, 5.9]))


class TestAsHandle(TestCase):
    """
    Test the asHandle function
    """

    def testOpenFile(self):
        """
        When an open file pointer is passed to asHandle, that same file
        pointer must be returned.
        """
        with patch.object(builtins, "open", mock_open()):
            fp = open("file")
            with asHandle(fp) as newfp:
                self.assertIs(fp, newfp)

    def testStr(self):
        """
        When a string filename is passed to asHandle, it must be possible to
        read the correct data from the fp that is returned.
        """
        mockOpener = mock_open(read_data="xxx")
        with patch.object(builtins, "open", mockOpener):
            with asHandle("file") as fp:
                self.assertEqual("xxx", fp.read())

    def testBZ2(self):
        """
        When a string '*.bz2' filename is passed to asHandle, it must be
        possible to read the correct data from the fp that is returned.
        """
        result = BytesIO(b"xxx")

        with patch.object(bz2, "BZ2File") as mockMethod:
            mockMethod.return_value = result
            with asHandle("file.bz2") as fp:
                self.assertEqual("xxx", fp.read())

    def testGzip(self):
        """
        When a string '*.gz' filename is passed to asHandle, it must be
        possible to read the correct data from the fp that is returned.
        """
        result = BytesIO(b"xxx")

        with patch.object(gzip, "GzipFile") as mockMethod:
            mockMethod.return_value = result
            with asHandle("file.gz") as fp:
                self.assertEqual("xxx", fp.read())


class TestParseRangeString(TestCase):
    """
    Check that the parseRangeString function works as expected.
    """

    def testEmptyString(self):
        """
        An empty string must produce an empty set of indices.
        """
        error = "^Illegal range ''. Ranges must single numbers or number-number\\.$"
        self.assertRaisesRegex(ValueError, error, parseRangeString, "")

    def testSingleNumber(self):
        """
        A single number must result in the expected set.
        """
        self.assertEqual({6}, parseRangeString("6"))

    def testSingleNumberSpaceBefore(self):
        """
        A single number preceeded by whitespace must result in the expected
        set.
        """
        self.assertEqual({6}, parseRangeString("  6"))

    def testSingleNumberSpaceAfter(self):
        """
        A single number followed by whitespace must result in the expected
        set.
        """
        self.assertEqual({6}, parseRangeString("6  "))

    def testSingleNumberSpaceBeforeAndAfter(self):
        """
        A single number preceeded and followed by whitespace must result in
        the expected set.
        """
        self.assertEqual({6}, parseRangeString(" 6  "))

    def testSingleRange(self):
        """
        A single range must result in the expected set.
        """
        self.assertEqual({6, 7, 8, 9, 10}, parseRangeString("6-10"))

    def testSingleRangeWithSpaceBeforeHyphen(self):
        """
        A single range with a space before the hyphen must result in the
        expected set.
        """
        self.assertEqual({6, 7, 8, 9, 10}, parseRangeString("6 -10"))

    def testSingleRangeWithSpaceAfterHyphen(self):
        """
        A single range with a space after the hyphen must result in the
        expected set.
        """
        self.assertEqual({6, 7, 8, 9, 10}, parseRangeString("6- 10"))

    def testSingleRangeWithSpaceBeforeAfterHyphen(self):
        """
        A single range with spaces before and after the hyphen must result in
        the expected set.
        """
        self.assertEqual({6, 7, 8, 9, 10}, parseRangeString("6 - 10"))

    def testTwoRanges(self):
        """
        Two ranges must result in the expected set.
        """
        self.assertEqual({6, 7, 8, 9, 10}, parseRangeString("6-8,9-10"))

    def testTwoOverlappingRanges(self):
        """
        Two overlapping ranges must result in the expected set.
        """
        self.assertEqual({6, 7, 8, 9, 10}, parseRangeString("6-9,7-10"))

    def testTwoRangesAndANumber(self):
        """
        Two ranges and a number must result in the expected set.
        """
        self.assertEqual({6, 7, 8, 10}, parseRangeString("6-8,10"))

    def testTwoRangesAndTwoNumbers(self):
        """
        Two ranges and two numbers must result in the expected set.
        """
        self.assertEqual({4, 6, 7, 8, 9, 10, 11, 12}, parseRangeString("6-8,9,10-12,4"))

    def testZeroConversion(self):
        """
        If we ask for zero conversion, the result must be as expected.
        """
        self.assertEqual(
            {3, 5, 6, 7, 8, 9, 10, 11},
            parseRangeString("6-8,9,10-12,4", convertToZeroBased=True),
        )

    def testAsListPreservesOrderAlreadyOrdered(self):
        """
        When asList is True, the returned value must be a list and the order
        of the numbers in the range must be preserved. Here the passed range
        is already in ascending order.
        """
        result = parseRangeString("6-9,7-10", asList=True)
        self.assertIsInstance(result, list)
        self.assertEqual([6, 7, 8, 9, 10], result)

    def testAsListPreservesOrderDisordered(self):
        """
        When asList is True, the returned value must be a list and the order
        of the numbers in the range must be preserved. Here the passed range
        is not in ascending order (and must be ordered identically in the result).
        """
        result = parseRangeString("7-10,1-3", asList=True)
        self.assertIsInstance(result, list)
        self.assertEqual([7, 8, 9, 10, 1, 2, 3], result)

    def testAsListIgnoresDuplicates(self):
        """
        When asList is True, the returned value must not include duplicates.
        """
        result = parseRangeString("7-10,1-3,8-12,2", asList=True)
        self.assertIsInstance(result, list)
        self.assertEqual([7, 8, 9, 10, 1, 2, 3, 11, 12], result)


class TestParseRangeExpression(TestCase):
    """
    Check that the parseRangeExpression function works as expected.
    """

    def testInvalidExpression(self):
        """
        An invalid string must raise a ValueError.
        """
        error = r"^\($"
        self.assertRaisesRegex(ValueError, error, parseRangeExpression, "(")
        error = r"^hey$"
        self.assertRaisesRegex(ValueError, error, parseRangeExpression, "hey")

    def testEmptyString(self):
        """
        An empty string must produce an empty set.
        """
        self.assertEqual(set(), parseRangeExpression(""))

    def testOneRange(self):
        """
        A simple 3-4 string must produce the expected set.
        """
        self.assertEqual({3, 4}, parseRangeExpression("3-4"))

    def testOneRangeZeroBased(self):
        """
        A simple 3-4 string must produce the expected set when
        convertToZeroBased is True.
        """
        self.assertEqual({2, 3}, parseRangeExpression("3-4", True))

    def testCommas(self):
        """
        A simple 3,4,5 string must produce the expected set.
        """
        self.assertEqual({3, 4, 5}, parseRangeExpression("3,4,5"))

    def testCommasAndRange(self):
        """
        A simple 3,4,5-7 string must produce the expected set.
        """
        self.assertEqual({3, 4, 5, 6, 7}, parseRangeExpression("3,4,5-7"))

    def testTwoRanges(self):
        """
        A simple 3-4,6-8 string must produce the expected set.
        """
        self.assertEqual({3, 4, 6, 7, 8}, parseRangeExpression("3-4,6-8"))

    def testTwoRangesWithSpace(self):
        """
        A simple 3-4, 6-8 string must produce the expected set.
        """
        self.assertEqual({3, 4, 6, 7, 8}, parseRangeExpression("3-4, 6-8"))

    def testUnion(self):
        """
        A union such as 3-4 | 6-8 must produce the expected set.
        """
        self.assertEqual({3, 4, 6, 7, 8}, parseRangeExpression("3-4 | 6-8"))

    def testIntersection(self):
        """
        An intersection such as 3-4 & 4-8 must produce the expected set.
        """
        self.assertEqual({4}, parseRangeExpression("3-4 & 4-8"))

    def testDifferenceNoSpaces(self):
        """
        A difference such as 6-10-7-8 must produce the expected set.
        """
        self.assertEqual({6, 9, 10}, parseRangeExpression("6-10-7-8"))

    def testDifferenceWithSpaces(self):
        """
        A difference such as 6-10 - 7-8 must produce the expected set.
        """
        self.assertEqual({6, 9, 10}, parseRangeExpression("6-10 - 7-8"))

    def testParens(self):
        """
        A difference with parentheses such as '(3-5 | 7-9) & 5-7' must produce
        the expected set.
        """
        self.assertEqual({5, 7}, parseRangeExpression("(3-5 | 7-9) & 5-7"))

    def testDoubleParens(self):
        """
        A difference with two parentheses such as '(3-5 | 7-9) & (5-7 | 9-11)'
        must produce the expected set.
        """
        self.assertEqual({5, 7, 9}, parseRangeExpression("(3-5 | 7-9) & (5-7 | 9-11)"))


class TestStringIO(TestCase):
    """
    Tests for our StringIO class.
    """

    def testInitiallyEmpty(self):
        """
        A StringIO instance must initially be empty.
        """
        self.assertEqual("", StringIO().getvalue())

    def testWriteRead(self):
        """
        It must be possible to write and read to/from a StringIO instance as
        normal.
        """
        s = StringIO()
        s.write("hey")
        self.assertEqual("hey", s.getvalue())

    def testInitializedRead(self):
        """
        It must be possible to read from a StringIO instance that is
        initialized on creation.
        """
        s = StringIO("hey")
        self.assertEqual("hey", s.getvalue())

    def testContextManager(self):
        """
        It must be possible to use a StringIO instance as a context manager.
        """
        with StringIO() as s:
            s.write("hey")
            self.assertEqual("hey", s.getvalue())


class TestBaseCountsToStr(TestCase):
    """
    Test the baseCountsToStr function.
    """

    def testSimple(self):
        """
        A simple example must work as expected.
        """
        counts = Counter()
        counts["A"] += 1
        counts["G"] += 2
        self.assertEqual("A:1 G:2", baseCountsToStr(counts))


class TestNucleotidesToStr(TestCase):
    """
    Test the nucleotidesToStr function.
    """

    def testSimple(self):
        """
        A simple example must work as expected.
        """
        counts1 = Counter()
        counts1["A"] += 1
        counts1["G"] += 2
        counts2 = Counter()
        counts2["C"] += 1
        counts2["T"] += 3
        self.assertEqual(
            "0: A:1 G:2\n7: C:1 T:3",
            nucleotidesToStr(
                {
                    0: counts1,
                    7: counts2,
                }
            ),
        )


class TestCountPrint(TestCase):
    """
    Test the countPrint function and the contained percentage function.
    """

    def testSimple(self):
        """
        A simple example must work as expected.
        """
        count = 2
        len1 = 10
        self.assertEqual("Count is: 2/10 (20.00%)", countPrint("Count is", count, len1))

    def testTwoSequences(self):
        """
        An example involving two different lengths must work as expected.
        """
        count = 2
        len1 = 10
        len2 = 8
        self.assertEqual(
            "Count is: 2/10 (20.00%) of sequence 1," " 2/8 (25.00%) of sequence 2",
            countPrint("Count is", count, len1, len2),
        )


class TestPct(TestCase):
    """
    Test the pct function.
    """

    def testZeroNumerator(self):
        """
        The pct function must produce the correct result if the numerator is
        zero.
        """
        self.assertEqual("0/10 (0.000%)", pct(0, 10))

    def testZeroDenominator(self):
        """
        The pct function must produce the correct result if the denominator is
        zero.
        """
        self.assertEqual("0/0 (0.000%)", pct(0, 0))

    def testOneHalf(self):
        """
        The pct function must produce the correct result if the numerator is
        one half of the denominator.
        """
        self.assertEqual("5/10 (50.000%)", pct(5, 10))

    def testOneSeventh(self):
        """
        The pct function must produce the correct result if the numerator is
        one seventh of the denominator.
        """
        self.assertEqual("2/14 (14.286%)", pct(2, 14))


class TestTake(TestCase):
    """
    Test the take function.
    """

    def testNLessThanOne(self):
        """
        The take function must raise an AssertionError if passed n < 1.
        list.
        """
        self.assertRaises(AssertionError, next, take([], 0))

    def testEmpty(self):
        """
        The take function must return an empty generator if passed an empty
        list.
        """
        self.assertEqual([], list(take([], 4)))

    def testOne(self):
        """
        The take function must return individual items if passed n=1.
        """
        self.assertEqual([[3], [4], [5]], list(take([3, 4, 5], 1)))

    def testNBiggerThanPassedList(self):
        """
        The take function must return all items if passed an n that is bigger
        than the passed list.
        """
        self.assertEqual([[3, 4, 5]], list(take([3, 4, 5], 100)))

    def testThree(self):
        """
        The take function must return three items at a time if passed n=3, and
        also return the extra two items.
        """
        self.assertEqual(
            [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10]], list(take(range(11), 3))
        )


class TestMatchOffset(TestCase):
    """
    Test the matchOffset function.
    """

    def testEmpty(self):
        """
        An empty query must match an empty reference at position 0.
        """
        self.assertEqual(0, matchOffset("", ""))

    def testEqualStrings(self):
        """
        An non-empty reference must match an identical non-empty query at
        position 0.
        """
        self.assertEqual(0, matchOffset("AA", "AA"))

    def testQueryPaddedLeftByOne(self):
        """
        A query that is padded on the left by one space must match at
        position 1.
        """
        self.assertEqual(1, matchOffset("AA", " A"))

    def testQueryPaddedLeftByTwo(self):
        """
        A query that is padded on the left by two spaces must match at
        position 2.
        """
        self.assertEqual(2, matchOffset("AAA", "  A"))

    def testQueryPaddedLeftByTwoReferencePaddedLeftByOne(self):
        """
        A query that is padded on the left by two spaces must match a reference
        that is padded on the left by one space at position 1.
        """
        self.assertEqual(1, matchOffset(" AA", "  A"))

    def testQueryPaddedLeftByFiveReferencePaddedLeftByOneWithGaps(self):
        """
        A query that is padded on the left by five spaces must correctly match
        a reference containing gaps that is padded on the left by one space.
        """
        self.assertEqual(2, matchOffset(" AA--G", "     A"))

    def testQueryPaddedLeftByFiveReferenceOneGapLeftWithGaps(self):
        """
        A query that is padded on the left by five spaces must correctly match
        a reference that starts with a gap and that contains gaps.
        """
        self.assertEqual(2, matchOffset("-AA--G", "     A"))

    def testOmicronPartialInsertionRead(self):
        """
        A query that overlaps part of an insertion in the reference must be
        handled correctly.
        """
        self.assertEqual(
            12,
            matchOffset(
                "TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
                "             AGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTT",
            ),
        )


class TestIntsToIntervals(TestCase):
    """
    Test the intsToIntervals function.
    """

    def testEmpty(self):
        "An empty set of numbers must give an empty result."
        self.assertEqual(intsToIntervals([]), [])

    def testOneNumber(self):
        "A set of one number must result in that number."
        self.assertEqual(intsToIntervals({2018}), [(2018, 2018)])

    def testTwoNonConsecutiveNumbers(self):
        """
        A set of two numbers that are not consecutive must result in those numbers.
        """
        self.assertEqual(intsToIntervals({2018, 2020}), [(2018, 2018), (2020, 2020)])

    def testTwoNonConsecutiveNumbersDecreasingOrder(self):
        """
        A set of two numbers that are not consecutive (and which decrease in the
        passed order) must result in those numbers.
        """
        self.assertEqual(intsToIntervals([2020, 2018]), [(2018, 2018), (2020, 2020)])

    def testTwoConsecutiveNumbers(self):
        """
        A set of two numbers that are consecutive must result in those numbers in a
        range.
        """
        self.assertEqual(intsToIntervals({2018, 2019}), [(2018, 2019)])

    def testThreeConsecutiveNumbers(self):
        """
        A set of three numbers that are consecutive must result in those numbers in a
        range.
        """
        self.assertEqual(intsToIntervals({2018, 2019, 2020}), [(2018, 2020)])

    def testTwoConsecutiveIntervals(self):
        """
        Two set of consecutive numbers must result in those numbers in two
        intervals.
        """
        self.assertEqual(
            intsToIntervals({2018, 2019, 2020, 2014, 2015}),
            [(2014, 2015), (2018, 2020)],
        )

    def testOneToAHundred(self):
        """
        1 to 100 must result in "1-100".
        """
        self.assertEqual(intsToIntervals(range(1, 101)), [(1, 100)])

    def testNegativeRange(self):
        """
        A range from one negative number to another must give the expected result.
        """
        self.assertEqual(intsToIntervals({-2018, -2019}), [(-2019, -2018)])

    def testPositiveAndNegativeIntervals(self):
        """
        Two intervals, one negative and one positive must give the expected result.
        """
        self.assertEqual(
            intsToIntervals({-2020, -2018, -2019, 2010, 2011, 2012}),
            [(-2020, -2018), (2010, 2012)],
        )


class TestIntsToStringIntervals(TestCase):
    """
    Test the intsToStringIntervals function.
    """

    def testEmpty(self):
        "An empty set of numbers must give an empty result."
        self.assertEqual(intsToStringIntervals([]), [])

    def testOneNumber(self):
        "A set of one number must result in that number."
        self.assertEqual(intsToStringIntervals({2018}), ["2018"])

    def testTwoNonConsecutiveNumbers(self):
        """
        A set of two numbers that are not consecutive must result in those numbers.
        """
        self.assertEqual(intsToStringIntervals({2018, 2020}), ["2018", "2020"])

    def testTwoNonConsecutiveNumbersDecreasingOrder(self):
        """
        A set of two numbers that are not consecutive (and which decrease in the
        passed order) must result in those numbers.
        """
        self.assertEqual(intsToStringIntervals([2020, 2018]), ["2018", "2020"])

    def testTwoConsecutiveNumbers(self):
        """
        A set of two numbers that are consecutive must result in those numbers in a
        range.
        """
        self.assertEqual(intsToStringIntervals({2018, 2019}), ["2018-2019"])

    def testTwoConsecutiveNumbersAlternateSeparator(self):
        """
        A set of two numbers that are consecutive must result in those numbers in a
        range, separated by the passed separator string.
        """
        self.assertEqual(
            intsToStringIntervals({2018, 2019}, sep=" -> "), ["2018 -> 2019"]
        )

    def testThreeConsecutiveNumbers(self):
        """
        A set of three numbers that are consecutive must result in those numbers in a
        range.
        """
        self.assertEqual(intsToStringIntervals({2018, 2019, 2020}), ["2018-2020"])

    def testTwoConsecutiveIntervals(self):
        """
        Two set of consecutive numbers must result in those numbers in two
        intervals.
        """
        self.assertEqual(
            intsToStringIntervals({2018, 2019, 2020, 2014, 2015}),
            ["2014-2015", "2018-2020"],
        )

    def testOneToAHundred(self):
        """
        1 to 100 must result in "1-100".
        """
        self.assertEqual(intsToStringIntervals(range(1, 101)), ["1-100"])

    def testNegativeRange(self):
        """
        A range from one negative number to another must give the expected result.
        """
        self.assertEqual(
            intsToStringIntervals({-2018, -2019}, sep=" -> "), ["-2019 -> -2018"]
        )

    def testPositiveAndNegativeIntervals(self):
        """
        Two intervals, one negative and one positive must give the expected result.
        """
        self.assertEqual(
            intsToStringIntervals({-2020, -2018, -2019, 2010, 2011, 2012}, sep=" -> "),
            ["-2020 -> -2018", "2010 -> 2012"],
        )

    def testTwoConsecutiveNumbersMinimized(self):
        """
        When minimize=True is passed, a set of two numbers that are consecutive must
        result in those numbers in a simplified range.
        """
        self.assertEqual(intsToStringIntervals({2018, 2019}, minimize=True), ["2018-9"])

    def testTwoConsecutiveNegativeNumbersMinimized(self):
        """
        When minimize=True is passed, a set of two negative numbers that are
        consecutive must result in those numbers in a simplified range.
        """
        self.assertEqual(
            intsToStringIntervals(
                {-2020, -2021, -2018, -2019}, minimize=True, sep=" -> "
            ),
            ["-2021 -> 18"],
        )

    def testThreeConsecutiveNumbersMinimized(self):
        """
        When minimize=True is passed, a set of three numbers that are consecutive must
        result in those numbers in a simplified range.
        """
        self.assertEqual(
            intsToStringIntervals({2018, 2019, 2020}, minimize=True), ["2018-20"]
        )

    def testNegativePositiveRangeIsNotMinimized(self):
        """
        When minimize=True is passed, a set of numbers that are consecutive and
        crossing zero must result in that number range not being minimized.
        """
        self.assertEqual(
            intsToStringIntervals(range(-20, 21), minimize=True, sep=" -> "),
            ["-20 -> 20"],
        )

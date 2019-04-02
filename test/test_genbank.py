from unittest import TestCase
import six

from dark.genbank import splitRange, splitBioPythonRange


class TestSplitRangeString(TestCase):
    """
    Test the splitRange function.
    """
    def testEmptyString(self):
        """
        The splitRange function must raise ValueError when passed an
        range empty string.
        """
        error = (r'^Could not parse GenBank range string ""\. Original '
                 r'parsing ValueError was "invalid literal for int\(\) with '
                 "base 10: ''\"\\.$")
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, '')

    def testUnparseableString(self):
        """
        The splitRange function must raise ValueError when passed an
        unparseable string.
        """
        error = (r'^Could not parse GenBank range string "x"\. Original '
                 r'parsing ValueError was "invalid literal for int\(\) with '
                 "base 10: 'x'\"\\.$")
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, 'x')

    def testJoinWithNoRange(self):
        """
        The splitRange function must raise ValueError when passed an
        'join()' string with no ranges.
        """
        error = (r'^Could not parse GenBank range string "join\(\)"\. '
                 r'join\(\) can only be used with multiple ranges\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, 'join()')

    def testComplementWithNoRange(self):
        """
        The splitRange function must raise ValueError when passed an
        'complement()' string.
        """
        error = (r'^Could not parse GenBank range string "complement\(\)"\. '
                 r'Original parsing ValueError was "invalid literal for '
                 r'int\(\) with '
                 "base 10: ''\"\\.$")
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, 'complement()')

    def testDecreasingOrder(self):
        """
        The splitRange function must raise ValueError when passed a string
        in which the start, stop values decrease.
        """
        error = (r'^Could not parse GenBank range string "7..4"\. '
                 r'Offset values \(7, 4\) cannot decrease\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, '7..4')

    def testTwoNakedRanges(self):
        """
        Two ranges with no join qualifier must result in a ValueError.
        """
        error = (r'^Could not parse GenBank range string "3..5,7..9"\. '
                 r'Multiple ranges must be wrapped in join\(\)\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, '3..5,7..9')

    def testOneJoinedRange(self):
        """
        A single range with a join() must result in a ValueError.
        """
        error = (r'^Could not parse GenBank range string "join\(3..5\)"\. '
                 r'join\(\) can only be used with multiple ranges\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, 'join(3..5)')

    def testTwoComplementedButNotJoinedRange(self):
        """
        Two complement ranges that are not joined must result in a ValueError.
        """
        error = (r'^Could not parse GenBank range string '
                 r'"complement\(3..5,7..9\)"\. '
                 r'Multiple ranges must be wrapped in join\(\)\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, 'complement(3..5,7..9)')

    def testJoinComplement(self):
        """
        A join(complement(...) argument must result in a ValueError.
        """
        error = (r'^Could not parse GenBank range string '
                 r'"join\(complement\(3..5,7..9\)\)"\. '
                 r'Original parsing ValueError was "invalid literal for '
                 r'int\(\) with '
                 "base 10: 'complement\\(3'\"\\.$")
        six.assertRaisesRegex(
            self, ValueError, error, splitRange, 'join(complement(3..5,7..9))')

    def testOneNakedRange(self):
        """
        A single range with no qualifiers must return the expected result.
        """
        self.assertEqual(
            (((2, 5),), False),
            splitRange('3..5'))

    def testTwoJoinedRange(self):
        """
        Two joined ranges must return the expected result.
        """
        self.assertEqual(
            (((2, 5), (6, 9)), False),
            splitRange('join(3..5,7..9)'))

    def testOneComplementRange(self):
        """
        A single range with a complement() must return the expected result.
        """
        self.assertEqual(
            (((2, 5),), True),
            splitRange('complement(3..5)'))

    def testTwoComplementRanges(self):
        """
        Two ranges with a complement() must return the expected result.
        """
        self.assertEqual(
            (((2, 5), (7, 20)), True),
            splitRange('complement(join(3..5,8..20))'))

    def testThreeComplementRanges(self):
        """
        Three ranges with a complement() must return the expected result.
        """
        self.assertEqual(
            (((2, 5), (7, 19), (29, 39)), True),
            splitRange('complement(join(3..5,8..19,30..39))'))


class TestSplitBioPythonRangeString(TestCase):
    """
    Test the splitBioPythonRange function.
    """
    def testEmptyString(self):
        """
        The splitBioPythonRange function must raise ValueError when passed an
        range empty string.
        """
        error = (r'^Could not parse GenBank range string ""\. '
                 r'Range "" does not end with \]\(\+\) or \]\(-\)\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, '')

    def testUnparseableString(self):
        """
        The splitBioPythonRange function must raise ValueError when passed an
        unparseable string.
        """
        error = (r'^Could not parse GenBank range string "x"\. '
                 r'Range "x" does not end with \]\(\+\) or \]\(-\)\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, 'x')

    def testNoComplementIndicator(self):
        """
        The splitBioPythonRange function must raise ValueError when passed an
        a range that does not end with ](+] or ](-).
        """
        error = (r'^Could not parse GenBank range string "\[33:40\]"\. '
                 r'Range "\[33:40\]" does not end with \]\(\+\) or \]\(-\)\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, '[33:40]')

    def testNoLeadingSquareBracket(self):
        """
        The splitBioPythonRange function must raise ValueError when passed an
        a range that does not start with [.
        """
        error = (r'^Could not parse GenBank range string "33:40\]\(-\)"\. '
                 r'Range "33:40\]\(-\)" does not start with "\["\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, '33:40](-)')

    def testJoinWithNoRange(self):
        """
        The splitBioPythonRange function must raise ValueError when passed an
        'join()' string with no ranges.
        """
        error = (r'^Could not parse GenBank range string "join{}"\. '
                 r'join{} can only be used with multiple ranges\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, 'join{}')

    def testDecreasingOrder(self):
        """
        The splitBioPythonRange function must raise ValueError when passed a
        string in which the start, stop values decrease.
        """
        error = (r'^Could not parse GenBank range string "\[7:4\]\(\+\)"\. '
                 r'Offset values \(7, 4\) cannot decrease\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, '[7:4](+)')

    def testTwoNakedRanges(self):
        """
        Two ranges with no join qualifier must result in a ValueError.
        """
        error = (r'^Could not parse GenBank range string '
                 r'"\[3:5\(\+\), 7:9\(-\)\]"\. '
                 r'Multiple ranges must be wrapped in join{}\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, '[3:5(+), 7:9(-)]')

    def testOneJoinedRange(self):
        """
        A single range with a join() must result in a ValueError.
        """
        error = (r'^Could not parse GenBank range string '
                 r'"join{\[3:5\]\(-\)}"\. '
                 r'join{} can only be used with multiple ranges\.$')
        six.assertRaisesRegex(
            self, ValueError, error, splitBioPythonRange, 'join{[3:5](-)}')

    def testOneNakedRange(self):
        """
        A single range with no qualifiers must return the expected result.
        """
        self.assertEqual(
            [(3, 5, False)],
            splitBioPythonRange('[3:5](+)'))

    def testTwoJoinedRange(self):
        """
        Two joined ranges must return the expected result.
        """
        self.assertEqual(
            [(3, 5, False), (7, 9, True)],
            splitBioPythonRange('join{[3:5](+), [7:9](-)}'))

    def testThreeJoinedRange(self):
        """
        Three joined ranges must return the expected result.
        """
        self.assertEqual(
            [(3, 5, False), (7, 9, True), (17, 19, True)],
            splitBioPythonRange('join{[3:5](+), [7:9](-), [17:19](-)}'))

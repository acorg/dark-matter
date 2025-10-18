from io import BytesIO
from unittest import TestCase

import pytest

from dark.analyze_reads import getPrefixAndSuffix, trimReads


class TestGetPrefixAndSuffix(TestCase):
    """
    Tests for getPrefixAndSuffix()
    """

    def testZeroInput(self):
        result = getPrefixAndSuffix(BytesIO())
        self.assertEqual(result, (0, 0))

    def testOneInput(self):
        seq = b">hey\nagtcagtcagtc"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (12, 12))

    def testTwoDifferentReads(self):
        seq = b">hey\nagtcagtcagtc\n>you\ntcctg"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (0, 0))

    def testSamePrefixTwoReads(self):
        seq = b">hey\nagtcagtcagtc\n>you\nagttcctg"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (3, 0))

    def testSameSuffixTwoReads(self):
        seq = b">hey\nagtcagtcagtc\n>you\ntcctggtc"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (0, 3))

    def testSamePrefixAndSuffixTwoReads(self):
        seq = b">hey\nagtcagtcagtc\n>you\nagttcctggtc"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (3, 3))

    def testSamePrefixDifferentSuffixThreeReads(self):
        seq = b">hey\nagtcagtcagtc\n>you\nagttcctggtc\n>how\nagtcggtat"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (3, 0))

    def testSamePrefixSameSuffixThreeReads(self):
        seq = b">hey\nagtccttagatcg\n>you\nagtcgaatcg\n>how\nagtaacctcg"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (3, 3))

    def testDifferentPrefixSameSuffixThreeReads(self):
        seq = b">hey\nagtccttagatcg\n>you\ncgaatcg\n>how\natgacctcg"
        result = getPrefixAndSuffix(BytesIO(seq))
        self.assertEqual(result, (0, 3))


#
# Tests for trimReads()
#
@pytest.mark.parametrize(
    "sequence,prefix,suffix,expected",
    (
        (b"", 0, 0, ""),
        (b"", 10, 10, ""),
        (b"agtccgatcg", 0, 0, "agtccgatcg"),
        (b"agtccgatcg", 3, 0, "ccgatcg"),
        (b"agtccgatcg", 0, 3, "agtccga"),
        (b"agtccgatcg", 3, 3, "ccga"),
        (b"agtccgatcg", 8, 4, ""),
        (b"agtccgatcg", 4, 8, ""),
    ),
)
def testTrimReads(sequence, prefix, suffix, expected):
    result = list(trimReads(prefix, suffix, BytesIO(b">id\n" + sequence + b"\n")))[0]
    assert result.sequence == expected

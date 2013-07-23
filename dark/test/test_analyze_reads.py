from unittest import TestCase
from dark.analyze_reads import getPrefixAndSuffix, trimReads
from cStringIO import StringIO
from Bio import SeqIO


class TestGetPrefixAndSuffix(TestCase):
    """
    Tests for analyze_reads.py
    """

    def testZeroInput(self):
        result = getPrefixAndSuffix(StringIO())
        self.assertEqual(result, (0, 0))

    def testOneInput(self):
        seq = '>hey\nagtcagtcagtc'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (12, 12))

    def testTwoDifferentReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\ntcctg'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (0, 0))

    def testSamePrefixTwoReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\nagttcctg'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (3, 0))

    def testSameSuffixTwoReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\ntcctggtc'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (0, 3))

    def testSamePrefixAndSuffixTwoReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\nagttcctggtc'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (3, 3))

    def testSamePrefixDifferentSuffixThreeReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\nagttcctggtc\n>how\nagtaatcggtac'
        result = getPrefixAndSuffix(SeqIO(seq))
        self.assertEqual(result, (3, 0))

    def testDifferentPrefixSameSuffixThreeReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\ntcctggtc\n>how\nagtaatcggtacgtc'
        result = getPrefixAndSuffix(SeqIO(seq))
        self.assertEqual(result, (0, 3))

    def testSamePrefixSameSuffixThreeReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\nagttcctggtc\n>how\nagtaatcggtacgtc'
        result = getPrefixAndSuffix(SeqIO(seq))
        self.assertEqual(result, (3, 3))


class TestTrimReads(TestCase):
    """
    Tests for trimReads
    """
    def testZeroPrefixZeroSuffix(self):
        seq = '>hey\nagtcagtcagtc'
        prefix = 0
        suffix = 0
        result = trimReads(prefix, suffix, SeqIO(seq))
        self.assertEqual(result, 'agtcagtcagtc')

    # tbc


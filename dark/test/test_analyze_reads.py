from unittest import TestCase
from dark.analyze_reads import getReads, _longestPrefixOfTwoSeqs
from cStringIO import StringIO
from Bio import SeqIO


class TestAnalyzeReads(TestCase):
    """
    Tests for analyze_reads.py
    """

    def testZeroInput(self):
        self.assertEqual(getReads("HCoV-EMC_0_reads.fasta"), [])

    def testOneInput(self):
        pass

    def testTwoDifferentReads(self):
        pass

    def testSamePrefixTwoReads(self):
        pass

    def testSameSuffixTwoReads(self):
        pass

    def testSamePrefixAndSuffixTwoReads(self):
        pass

    def testSamePrefixDifferentSuffixThreeReads(self):
        pass

    def testDifferentPrefixSameSuffixThreeReads(self):
        pass

    def testReverse(self):
        pass



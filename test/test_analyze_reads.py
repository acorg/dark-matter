from unittest import TestCase
from dark.analyze_reads import getPrefixAndSuffix, trimReads
from Bio import SeqIO

from dark.utils import StringIO


class TestGetPrefixAndSuffix(TestCase):
    """
    Tests for getPrefixAndSuffix()
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
        seq = '>hey\nagtcagtcagtc\n>you\nagttcctggtc\n>how\nagtcggtat'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (3, 0))

    def testSamePrefixSameSuffixThreeReads(self):
        seq = '>hey\nagtccttagatcg\n>you\nagtcgaatcg\n>how\nagtaacctcg'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (3, 3))

    def testDifferentPrefixSameSuffixThreeReads(self):
        seq = '>hey\nagtccttagatcg\n>you\ncgaatcg\n>how\natgacctcg'
        result = getPrefixAndSuffix(StringIO(seq))
        self.assertEqual(result, (0, 3))


class TestTrimReads(TestCase):
    """
    Tests for trimReads()
    """
    def testZeroInput(self):
        result = trimReads(0, 0, StringIO())
        test_result = list(SeqIO.parse(StringIO(), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

    def testOneInput(self):
        result = trimReads(10, 10, StringIO())
        test_result = list(SeqIO.parse(StringIO(), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

    def testZeroPrefixZeroSuffix(self):
        seq = '>hey\nagtccgatcg'
        result = trimReads(0, 0, StringIO(seq))
        trimmed_seq = '>hey\nagtccgatcg'
        test_result = list(SeqIO.parse(StringIO(trimmed_seq), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

    def testPrefixZeroSuffix(self):
        seq = '>hey\nagtccgatcg'
        trimmed_seq = '>hey\nccgatcg'
        result = trimReads(3, 0, StringIO(seq))
        test_result = list(SeqIO.parse(StringIO(trimmed_seq), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

    def testZeroPrefixSuffix(self):
        seq = '>hey\nagtccgatcg'
        trimmed_seq = '>hey\nagtccga'
        result = trimReads(0, 3, StringIO(seq))
        test_result = list(SeqIO.parse(StringIO(trimmed_seq), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

    def testPrefixSuffix(self):
        seq = '>hey\nagtccgatcg'
        trimmed_seq = '>hey\nccga'
        result = trimReads(3, 3, StringIO(seq))
        test_result = list(SeqIO.parse(StringIO(trimmed_seq), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

    def testLongPrefixSuffix(self):
        seq = '>hey\nagtccgatcg'
        trimmed_seq = '>hey\n'
        result = trimReads(8, 4, StringIO(seq))
        test_result = list(SeqIO.parse(StringIO(trimmed_seq), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

    def testPrefixLongSuffix(self):
        seq = '>hey\nagtccgatcg'
        trimmed_seq = '>hey\n'
        result = trimReads(4, 8, StringIO(seq))
        test_result = list(SeqIO.parse(StringIO(trimmed_seq), 'fasta'))
        self.assertEqual(list(map(str, tuple(result))),
                         list(map(str, tuple(test_result))))

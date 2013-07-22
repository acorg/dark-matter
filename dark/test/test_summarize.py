from unittest import TestCase
from cStringIO import StringIO
from collections import defaultdict
from Bio import SeqIO
from dark.summarize import summarize_reads


class TestSummarizeReads(TestCase):
    """
    Test the summarize_reads function.
    """

    def testReadNumberEmptyInput(self):
        result = summarize_reads(StringIO())
        self.assertEqual(result['read_number'], 0)

    def testReadNumberOneSequenceCount(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['read_number'], 1)

    def testReadNumberTwoSequencesCount(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['read_number'], 2)

    def testTotalLengthEmptyInput(self):
        result = summarize_reads(StringIO())
        self.assertEqual(result['total_length'], 0)

    def testTotalLengthOneString(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['total_length'], 12)

    def testTotalLengthTwoStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['total_length'], 17)

    def testBaseCountsEmptyImput(self):
        result = summarize_reads(StringIO())
        self.assertEqual(result['base_counts'], {})

    def testBaseCountsOneRead(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['base_counts'],
                         {'a': 3, 'c': 3, 't': 3, 'g': 3})

    def testBaseCountsTwoReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['base_counts'],
                         {'a': 4, 'c': 5, 't': 4, 'g': 4})

    def testMaxLengthListEmptyInput(self):
        result = summarize_reads(StringIO())
        self.assertEqual(result['max_length'], 0)

    def testMaxLengthListTwoStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['max_length'], 12)

    def testMinLengthListEmptyInput(self):
        result = summarize_reads(StringIO())
        self.assertEqual(result['min_length'], 0)

    def testMinLengthListTwoStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['min_length'], 5)

    def testSequenceListEmptyInput(self):
        result = summarize_reads(StringIO(), returnSequences=True)
        self.assertEqual(result['sequences'], [])

    def testSequenceListOneString(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarize_reads(StringIO(seq), returnSequences=True)
        test_result = list(SeqIO.parse(StringIO(seq), 'fasta'))
        self.assertEqual(result['sequences'], test_result)

    def testSequenceListTwoStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarize_reads(StringIO(seq), returnSequences=True)
        test_result = list(SeqIO.parse(StringIO(seq), 'fasta'))
        self.assertEqual(result['sequences'], test_result)

    def testMedianEmptyInput(self):
        result = summarize_reads(StringIO())
        self.assertEqual(result['median_length'], 0)

    def testMedianOneString(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['median_length'], 12)

    def testMedianThreeStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['median_length'], 7)

    def testMedianFourStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
        atggctattgaactgtatct'
        result = summarize_reads(StringIO(seq))
        self.assertEqual(result['median_length'], 9.5)

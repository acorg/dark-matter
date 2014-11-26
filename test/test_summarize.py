from unittest import TestCase
from cStringIO import StringIO

from dark.summarize import summarizeReads, summarizePosition
from dark.reads import Read, Reads


class TestSummarizeReads(TestCase):
    """
    Test the summarizeReads function.
    """

    def testReadNumberEmptyInput(self):
        result = summarizeReads(StringIO(), 'fasta')
        self.assertEqual(result['read_number'], 0)

    def testReadNumberOneSequenceCount(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['read_number'], 1)

    def testReadNumberTwoSequencesCount(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['read_number'], 2)

    def testTotalLengthEmptyInput(self):
        result = summarizeReads(StringIO(), 'fasta')
        self.assertEqual(result['total_length'], 0)

    def testTotalLengthOneString(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['total_length'], 12)

    def testTotalLengthTwoStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['total_length'], 17)

    def testBaseCountsEmptyImput(self):
        result = summarizeReads(StringIO(), 'fasta')
        self.assertEqual(result['base_counts'], {})

    def testBaseCountsOneRead(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['base_counts'],
                         {'a': 3, 'c': 3, 't': 3, 'g': 3})

    def testBaseCountsTwoReads(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['base_counts'],
                         {'a': 4, 'c': 5, 't': 4, 'g': 4})

    def testMaxLengthListEmptyInput(self):
        result = summarizeReads(StringIO(), 'fasta')
        self.assertEqual(result['max_length'], 0)

    def testMaxLengthListTwoStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['max_length'], 12)

    def testMinLengthListEmptyInput(self):
        result = summarizeReads(StringIO(), 'fasta')
        self.assertEqual(result['min_length'], 0)

    def testMinLengthListTwoStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['min_length'], 5)

    def testMedianEmptyInput(self):
        result = summarizeReads(StringIO(), 'fasta')
        self.assertEqual(result['median_length'], 0)

    def testMedianOneString(self):
        seq = '>hey\nagtcagtcagtc'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['median_length'], 12)

    def testMedianThreeStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['median_length'], 7)

    def testMedianFourStrings(self):
        seq = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
        atggctattgaactgtatct'
        result = summarizeReads(StringIO(seq), 'fasta')
        self.assertEqual(result['median_length'], 9.5)


class TestSummarizePosition(TestCase):
    """
    Tests for the summarizePosition function.
    """
    def testOffset(self):
        """
        The function must return the base at the right offset.
        """
        reads = Reads()
        reads.add(Read('hey', 'aataaa'))
        reads.add(Read('you', 'aata'))
        reads.add(Read('how', 'aataaaaaa'))
        result = summarizePosition(reads, 3)
        self.assertEqual(result['countAtPosition'], {'t': 3})

    def testNumberOfSequences(self):
        """
        Must return the right number of sequences.
        """
        reads = Reads()
        reads.add(Read('hey', 'agtcagtcagtc'))
        reads.add(Read('you', 'acctg'))
        reads.add(Read('how', 'atg'))
        result = summarizePosition(reads, 10)
        self.assertEqual(result['includedSequences'], 1)
        self.assertEqual(result['allSequences'], 3)

    def testFrequencies(self):
        """
        Must return the right frequencies.
        """
        reads = Reads()
        reads.add(Read('hey', 'aaaaaa'))
        reads.add(Read('you', 'aaca'))
        reads.add(Read('how', 'aataaaaaa'))
        result = summarizePosition(reads, 3)
        self.assertEqual(result['countAtPosition'], {'a': 1, 'c': 1, 't': 1})

    def testFrequenciesNonePresent(self):
        """
        Must return none, if no reads are present.
        """
        reads = Reads()
        result = summarizePosition(reads, 3)
        self.assertEqual(result['countAtPosition'], {})

    def testNumberOfSequencesNonePresent(self):
        """
        Must return none, if no reads are present.
        """
        reads = Reads()
        result = summarizePosition(reads, 3)
        self.assertEqual(result['includedSequences'], 0)

    def testPositionBiggerThanSequenceLength(self):
        """
        Must not count sequence if too short.
        """
        reads = Reads()
        reads.add(Read('hey', 'aaaaaa'))
        reads.add(Read('you', 'aaca'))
        reads.add(Read('how', 'aat'))
        result = summarizePosition(reads, 6)
        self.assertEqual(result['countAtPosition'], {'a': 1})

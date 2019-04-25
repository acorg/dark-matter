from unittest import TestCase
from six import assertRaisesRegex

from dark.reads import DNARead
from dark.summarize import summarizeReads, sequenceCategoryLengths
from dark.utils import StringIO


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


class TestSequenceCategoryLengths(TestCase):
    """
    Test the sequenceCategoryLengths function.
    """
    def testInvalidMinLength(self):
        """
        If a minLength value less than 1 is passed, a ValueError must be
        raised.
        """
        read = DNARead('id', '')
        error = '^minLength must be at least 1$'
        assertRaisesRegex(self, ValueError, error, sequenceCategoryLengths,
                          read, {}, minLength=0)

    def testEmpty(self):
        """
        An empty sequence should result in an empty category summary.
        """
        read = DNARead('id', '')
        self.assertEqual([], sequenceCategoryLengths(read, {}))

    def testOneCategoryPerBase(self):
        """
        If each base is in its own category, the summary must be correct.
        """
        read = DNARead('id', 'ACGT')
        categories = {
            'A': 0,
            'C': 1,
            'G': 2,
            'T': 3,
        }
        self.assertEqual([(0, 1), (1, 1), (2, 1), (3, 1)],
                         sequenceCategoryLengths(read, categories))

    def testRepeatedCategory(self):
        """
        If categories are repeated in a sequence, the summary must have the
        correct length for the categories.
        """
        read = DNARead('id', 'ACCGGTTT')
        categories = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('a', 1), ('c', 2), ('g', 2), ('t', 3)],
                         sequenceCategoryLengths(read, categories))

    def testUnknownCategory(self):
        """
        If a base has no category, the summary must have C{None} as the
        category for those bases.
        """
        read = DNARead('id', 'ACCGGTTT')
        categories = {
            'A': 'a',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('a', 1), (None, 2), ('g', 2), ('t', 3)],
                         sequenceCategoryLengths(read, categories))

    def testUnknownCategoryWithDefault(self):
        """
        If a base has no category, the summary must have the passed default
        category as the category for those bases.
        """
        read = DNARead('id', 'ACCGGTTT')
        categories = {
            'A': 'a',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('a', 1), ('xxx', 2), ('g', 2), ('t', 3)],
                         sequenceCategoryLengths(read, categories, 'xxx'))

    def testSuppressAtStart(self):
        """
        If a region at the start of the sequence is shorter than the passed
        minimum length, the result should suppress the catgeory information.
        """
        read = DNARead('id', 'ACCGGTTT')
        categories = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('...', 1), ('c', 2), ('g', 2), ('t', 3)],
                         sequenceCategoryLengths(read, categories,
                                                 minLength=2))

    def testSuppressTwoAtStart(self):
        """
        If 2 regions at the start of the sequence are shorter than the passed
        minimum length, the result should suppress the catgeory information
        and the length of the suppressed region must be the sum of the lengths
        of the regions.
        """
        read = DNARead('id', 'AGCCGGTTT')
        categories = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('...', 2), ('c', 2), ('g', 2), ('t', 3)],
                         sequenceCategoryLengths(read, categories,
                                                 minLength=2))

    def testSuppressAtEnd(self):
        """
        If a region at the end of the sequence is shorter than the passed
        minimum length, the result should suppress the catgeory information.
        """
        read = DNARead('id', 'CCGGTTTA')
        categories = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('c', 2), ('g', 2), ('t', 3), ('...', 1)],
                         sequenceCategoryLengths(read, categories,
                                                 minLength=2))

    def testSuppressTwoAtEnd(self):
        """
        If 2 regions at the end of the sequence are shorter than the passed
        minimum length, the result should suppress the catgeory information
        and the length of the suppressed region must be the sum of the lengths
        of the regions.
        """
        read = DNARead('id', 'CCGGTTTAC')
        categories = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('c', 2), ('g', 2), ('t', 3), ('...', 2)],
                         sequenceCategoryLengths(read, categories,
                                                 minLength=2))

    def testSuppressWithNonDefaultSuppresscategory(self):
        """
        If a region of the sequence is shorter than the passed minimum length,
        the result should suppress the catgeory information and the suppress
        category returned must be the one that is passed.
        """
        read = DNARead('id', 'ACCGGTTT')
        categories = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('s', 1), ('c', 2), ('g', 2), ('t', 3)],
                         sequenceCategoryLengths(
                             read, categories, minLength=2,
                             suppressedCategory='s'))

    def testAllSuppressed(self):
        """
        If all regions of the sequence are shorter than the passed
        minimum length, the result should suppress the catgeory information
        and the suppressed region length must be the sum of the region lengths.
        """
        read = DNARead('id', 'ACCGGGTTT')
        categories = {
            'A': 'a',
            'C': 'c',
            'G': 'g',
            'T': 't',
        }
        self.assertEqual([('...', 9)],
                         sequenceCategoryLengths(read, categories,
                                                 minLength=5))

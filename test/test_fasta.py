from dark.fasta import dedupFasta, dePrefixAndSuffixFasta
from cStringIO import StringIO
from unittest import TestCase
from Bio import SeqIO


class FastaDeDup(TestCase):
    """Tests for de-duping FASTA sequence lists."""

    def testEmpty(self):
        """An empty FASTA list gets de-duped to an empty list."""
        self.assertEqual(list(dedupFasta([])), [])

    def testLengthOne(self):
        """A FASTA list with just one item gets de-duped to the same one item."""
        seq = '>hey\nagtcagtcagtc'
        s1 = SeqIO.read(StringIO(seq), 'fasta')
        self.assertEqual(list(dedupFasta([s1])), [s1])

    def testRemovalOfIdenticalSequences(self):
        """A list with 2 copies of the same seq is de-duped to have 1 copy."""
        seq = '>hey\nagtcagtcagtc'
        s1 = SeqIO.read(StringIO(seq), 'fasta')
        s2 = SeqIO.read(StringIO(seq), 'fasta')
        self.assertEqual(list(dedupFasta([s1, s2])), [s1])


class Unused(TestCase):

    def testEmpty(self):
        """An empty FASTA list gets de-duped to an empty list."""
        self.assertEqual(list(dePrefixAndSuffixFasta([])), [])

    def testLengthOne(self):
        """A FASTA list with just one item gets de-duped to the same one item."""
        seq = '>hey\nagtcagtcagtc'
        s1 = SeqIO.read(StringIO(seq), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1])), [s1])

    def testRemovalOfIdenticalSequences(self):
        """A list with 2 copies of the same seq is de-duped to have 1 copy."""
        seq = '>hey\nagtcagtcagtc'
        s1 = SeqIO.read(StringIO(seq), 'fasta')
        s2 = SeqIO.read(StringIO(seq), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2])), [s1])

    def testRemovalOfPrefix(self):
        """A sequence that is a prefix of another is removed."""
        s1 = SeqIO.read(StringIO('>s1\nagtcagtcagtc'), 'fasta')
        s2 = SeqIO.read(StringIO('>s2\nagtcag'), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2])), [s1])

    def testRemovalOfSuffix(self):
        """A sequence that is a suffix of another is removed."""
        s1 = SeqIO.read(StringIO('>s1\nagtcagtcagtc'), 'fasta')
        s2 = SeqIO.read(StringIO('>s2\ncagtc'), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2])), [s1])

    def testRemovalOfPrefixSuffixAndDuplicate(self):
        """Prefixes, suffixes, and duplicates should collectively all be removed."""
        s1 = SeqIO.read(StringIO('>s1\nagtcagtcagtc'), 'fasta')
        s2 = SeqIO.read(StringIO('>s2\nagtcagtcagtc'), 'fasta')
        s3 = SeqIO.read(StringIO('>s3\nagtcagt'), 'fasta')
        s4 = SeqIO.read(StringIO('>s4\ntcagtc'), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2, s3, s4])), [s1])

    def testOrderIndependent(self):
        """A sequence that is a prefix of another is removed when it appears first."""
        s1 = SeqIO.read(StringIO('>s1\nagtcag'), 'fasta')
        s2 = SeqIO.read(StringIO('>s2\nagtcagtcagtc'), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2])), [s2])

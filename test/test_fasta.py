from dark.reads import Read
from dark.fasta import (dedupFasta, dePrefixAndSuffixFasta, fastaSubtract,
                        FastaReads)

from cStringIO import StringIO
from unittest import TestCase
from Bio import SeqIO
from mock import patch
from mocking import mockOpen


class FastaDeDup(TestCase):
    """Tests for de-duping FASTA sequence lists."""

    def testEmpty(self):
        """An empty FASTA list gets de-duped to an empty list."""
        self.assertEqual(list(dedupFasta([])), [])

    def testLengthOne(self):
        """
        A FASTA list with just one item gets de-duped to the same one item.
        """
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
        """
        A FASTA list with just one item gets de-duped to the same one item.
        """
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
        """
        Prefixes, suffixes, and duplicates should collectively all be removed.
        """
        s1 = SeqIO.read(StringIO('>s1\nagtcagtcagtc'), 'fasta')
        s2 = SeqIO.read(StringIO('>s2\nagtcagtcagtc'), 'fasta')
        s3 = SeqIO.read(StringIO('>s3\nagtcagt'), 'fasta')
        s4 = SeqIO.read(StringIO('>s4\ntcagtc'), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2, s3, s4])), [s1])

    def testOrderIndependent(self):
        """
        A sequence that is a prefix of another is removed when it appears
        first.
        """
        s1 = SeqIO.read(StringIO('>s1\nagtcag'), 'fasta')
        s2 = SeqIO.read(StringIO('>s2\nagtcagtcagtc'), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2])), [s2])


class TestFastaSubtract(TestCase):
    """
    Test the fastaSubtract() function.
    """
    def testZeroFiles(self):
        self.assertRaises(IndexError, fastaSubtract, [])

    def testOneFile(self):
        """
        When just one file is passed we should get a result that has as many
        reads as was in the single input file.
        """
        fasta1 = '\n'.join([
            '>one',
            'agtcagtcagtc',
            '>two',
            'acctg',
            '>three',
            'atgggtc',
            '>four',
            'atggctattgaactgtatct',
            ])

        result = list(fastaSubtract([StringIO(fasta1)]))
        self.assertEqual(len(result), 4)

    def testSubtractEverything(self):
        """
        When two input files have the same reads, subtraction must result in an
        empty (no reads) output.
        """
        fasta1 = '\n'.join([
            '>one',
            'agtcagtcagtc',
            '>two',
            'acctg',
            '>three',
            'atgggtc',
            '>four',
            'atggctattgaactgtatct',
            ])

        result = list(fastaSubtract([StringIO(fasta1), StringIO(fasta1)]))
        self.assertEqual([], result)

    def testSubtractFromNothing(self):
        """
        When the first file is empty, the result shoud be too.
        """
        fasta1 = ''
        fasta2 = '\n'.join([
            '>five',
            'agtcagtcagtc',
            '>six',
            'acctg',
            ])

        result = list(fastaSubtract([StringIO(fasta1), StringIO(fasta2)]))
        self.assertEqual([], result)

    def testSubtractNothing(self):
        """
        When two input files have no overlap, subtraction must result in the
        same reads as are in the first input.
        """
        fasta1 = '\n'.join([
            '>one',
            'agtcagtcagtc',
            '>two',
            'acctg',
            '>three',
            'atgggtc',
            '>four',
            'atggctattgaactgtatct',
            ])
        fasta2 = '\n'.join([
            '>five',
            'agtcagtcagtc',
            '>six',
            'acctg',
            ])

        result = list(fastaSubtract([StringIO(fasta1), StringIO(fasta2)]))
        self.assertEqual(['four', 'one', 'three', 'two'],
                         sorted([seq.id for seq in result]))

    def testThreeFiles(self):
        """
        Subtraction of three files must work correctly.
        """
        fasta1 = '\n'.join([
            '>one',
            'agtcagtcagtc',
            '>two',
            'acctg',
            '>three',
            'atgggtc',
            '>four',
            'atggctattgaactgtatct',
            ])
        fasta2 = '\n'.join([
            '>one',
            'agtcagtcagtc',
            ])
        fasta3 = '\n'.join([
            '>two',
            'acctg',
            '>three',
            'atgggtc',
            ])

        result = list(fastaSubtract([StringIO(fasta1),
                                     StringIO(fasta2),
                                     StringIO(fasta3)]))
        self.assertEqual(len(result), 1)
        self.assertEqual(str(result[0].seq), 'atggctattgaactgtatct')
        self.assertEqual(str(result[0].id), 'four')

    def testSequencesAreChecked(self):
        """
        If a two reads with the same id do not have the same sequence,
        an assertion error must be raised.
        """
        fasta1 = '\n'.join([
            '>one',
            'ag',
            ])
        fasta2 = '\n'.join([
            '>one',
            'at',
            ])

        self.assertRaises(AssertionError, fastaSubtract,
                          [StringIO(fasta1), StringIO(fasta2)])


class TestFastaReads(TestCase):
    """Tests for the L{dark.fasta.FastaReads} class."""

    def testEmpty(self):
        """An empty FASTA file results in an empty iterator."""
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            reads = FastaReads('filename.fasta')
            self.assertEqual([], list(reads))

    def testOneRead(self):
        """
        A FASTA file with one read must be read properly.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch('__builtin__.open', mockOpener, create=True):
            reads = list(FastaReads('filename.fasta'))
            self.assertEqual([Read('id1', 'ACGT')], reads)

    def testNoQuality(self):
        """
        A FASTA file read must not have any quality information.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch('__builtin__.open', mockOpener, create=True):
            reads = list(FastaReads('filename.fasta'))
            self.assertEqual(None, reads[0].quality)

    def testTwoReads(self):
        """
        A FASTA file with two reads must be read properly and its
        sequences must be returned in the correct order.
        """
        data = '\n'.join(['>id1', 'ACGT', '>id2', 'TGCA'])
        mockOpener = mockOpen(read_data=data)
        with patch('__builtin__.open', mockOpener, create=True):
            reads = list(FastaReads('filename.fasta'))
            self.assertEqual(2, len(reads))
            self.assertEqual([Read('id1', 'ACGT'), Read('id2', 'TGCA')], reads)

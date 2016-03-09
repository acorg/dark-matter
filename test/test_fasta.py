import six
from six.moves import builtins
from six import StringIO
from unittest import TestCase
from Bio import SeqIO

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import mockOpen

from dark.reads import Read, AARead, DNARead, RNARead, Reads
from dark.fasta import (dedupFasta, dePrefixAndSuffixFasta, fastaSubtract,
                        FastaReads, combineReads)


class FastaDeDup(TestCase):
    """
    Tests for de-duping FASTA sequence lists.
    """

    def testEmpty(self):
        """
        An empty FASTA list gets de-duped to an empty list.
        """
        self.assertEqual(list(dedupFasta([])), [])

    def testLengthOne(self):
        """
        A FASTA list with just one item gets de-duped to the same one item.
        """
        reads = Reads()
        reads.add(Read('id', 'GGG'))
        self.assertEqual(list(dedupFasta(reads)), [Read('id', 'GGG')])

    def testRemovalOfIdenticalSequences(self):
        """
        A list with 2 copies of the same seq is de-duped to have 1 copy.
        """
        reads = Reads()
        reads.add(Read('id', 'GGG'))
        reads.add(Read('id', 'GGG'))
        self.assertEqual(list(dedupFasta(reads)), [Read('id', 'GGG')])

    def testRemovalOfIdenticalSequencesWithDifferingIds(self):
        """
        A list with 2 copies of the same seq is de-duped to have 1 copy,
        including when the read ids differ.
        """
        reads = Reads()
        reads.add(Read('id1', 'GGG'))
        reads.add(Read('id2', 'GGG'))
        self.assertEqual(list(dedupFasta(reads)), [Read('id1', 'GGG')])


class Unused(TestCase):

    def testEmpty(self):
        """
        An empty FASTA list gets de-duped to an empty list.
        """
        self.assertEqual(list(dePrefixAndSuffixFasta([])), [])

    def testLengthOne(self):
        """
        A FASTA list with just one item gets de-duped to the same one item.
        """
        seq = '>hey\nagtcagtcagtc'
        s1 = SeqIO.read(StringIO(seq), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1])), [s1])

    def testRemovalOfIdenticalSequences(self):
        """
        A list with 2 copies of the same seq is de-duped to have 1 copy.
        """
        seq = '>hey\nagtcagtcagtc'
        s1 = SeqIO.read(StringIO(seq), 'fasta')
        s2 = SeqIO.read(StringIO(seq), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2])), [s1])

    def testRemovalOfPrefix(self):
        """
        A sequence that is a prefix of another is removed.
        """
        s1 = SeqIO.read(StringIO('>s1\nagtcagtcagtc'), 'fasta')
        s2 = SeqIO.read(StringIO('>s2\nagtcag'), 'fasta')
        self.assertEqual(list(dePrefixAndSuffixFasta([s1, s2])), [s1])

    def testRemovalOfSuffix(self):
        """
        A sequence that is a suffix of another is removed.
        """
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
    """
    Tests for the L{dark.fasta.FastaReads} class.
    """

    def testEmpty(self):
        """
        An empty FASTA file results in an empty iterator.
        """
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            reads = FastaReads('filename.fasta')
            self.assertEqual([], list(reads))

    def testOneRead(self):
        """
        A FASTA file with one read must be read properly.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta'))
            self.assertEqual([Read('id1', 'ACGT')], reads)

    def testNoQuality(self):
        """
        A FASTA file read must not have any quality information.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta'))
            self.assertEqual(None, reads[0].quality)

    def testTwoReads(self):
        """
        A FASTA file with two reads must be read properly and its
        sequences must be returned in the correct order.
        """
        data = '\n'.join(['>id1', 'ACGT', '>id2', 'TGCA'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta'))
            self.assertEqual(2, len(reads))
            self.assertEqual([Read('id1', 'ACGT'), Read('id2', 'TGCA')], reads)

    def testTypeDefaultsToDNA(self):
        """
        A FASTA file whose type is not specified must result in reads that
        are instances of DNARead.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta'))
            self.assertTrue(isinstance(reads[0], DNARead))

    def testTypeAA(self):
        """
        A FASTA file whose read class is AARead must result in reads that
        are instances of AARead.
        """
        data = '\n'.join(['>id1', 'ACGST'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta', AARead))
            self.assertTrue(isinstance(reads[0], AARead))

    def testTypeDNA(self):
        """
        A FASTA file whose read class is DNARead must result in reads that
        are instances of DNARead.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta', DNARead))
            self.assertTrue(isinstance(reads[0], DNARead))

    def testTypeRNA(self):
        """
        A FASTA file whose read class is RNARead must result in reads that
        are instances of RNARead.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta', RNARead))
            self.assertTrue(isinstance(reads[0], RNARead))

    def testAlphabetIsCheckedAndRaisesValueErrorOnFirstRead(self):
        """
        The default behavior of a FastaReads instance is to check to ensure
        its sequences have the correct alphabet and to raise ValueError if not.
        A non-alphabetic character in the first read must be detected.
        """
        data = StringIO('\n'.join([
            '>one',
            'at-at',
        ]))
        error = ("^Read alphabet \('-AT'\) is not a subset of expected "
                 "alphabet \('ACDEFGHIKLMNPQRSTVWY'\) for read class "
                 "AARead\.$")
        six.assertRaisesRegex(self, ValueError, error, list,
                              FastaReads(data, AARead))

    def testAlphabetIsCheckedAndRaisesValueErrorOnSecondRead(self):
        """
        The default behavior of a FastaReads instance is to check to ensure
        its sequences have the correct alphabet and to raise ValueError if not.
        A non-alphabetic character in the second read must be detected.
        """
        data = StringIO('\n'.join([
            '>one',
            'atat',
            '>two',
            'at-at',
        ]))
        error = ("^Read alphabet \('-AT'\) is not a subset of expected "
                 "alphabet \('ACDEFGHIKLMNPQRSTVWY'\) for read class "
                 "AARead\.$")
        six.assertRaisesRegex(self, ValueError, error, list,
                              FastaReads(data, AARead))

    def testDisableAlphabetChecking(self):
        """
        It must be possible to have a FastaReads instance not do alphabet
        checking, if requested (by passing checkAlphabet=0).
        """
        data = StringIO('\n'.join([
            '>one',
            'at-at',
        ]))
        self.assertEqual(1, len(list(FastaReads(data, AARead,
                                                checkAlphabet=0))))

    def testOnlyCheckSomeAlphabets(self):
        """
        It must be possible to have the alphabets of only a certain number of
        reads checked. A non-alphabetic character in a later read must not
        stop that read from being processed.
        """
        data = StringIO('\n'.join([
            '>one',
            'atat',
            '>two',
            'at-at',
        ]))
        reads = list(FastaReads(data, AARead, checkAlphabet=1))
        self.assertEqual(2, len(reads))
        self.assertEqual('at-at', reads[1].sequence)

    def testConvertLowerToUpperCaseIfSpecifiedAARead(self):
        """
        A read needs to be converted from lower to upper case if specified.
        """
        data = '\n'.join(['>id1', 'actgs'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta', readClass=AARead,
                         upperCase=True))
            self.assertEqual([AARead('id1', 'ACTGS')], reads)

    def testConvertLowerToUpperCaseIfSpecifiedDNARead(self):
        """
        A read needs to be converted from lower to upper case if specified.
        """
        data = '\n'.join(['>id1', 'actg'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta', upperCase=True))
            self.assertEqual([AARead('id1', 'ACTG')], reads)

    def testDontConvertLowerToUpperCaseIfNotSpecified(self):
        """
        A read must not be converted from lower to upper case if not specified.
        """
        data = '\n'.join(['>id1', 'actgs'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastaReads('filename.fasta', readClass=AARead))
            self.assertEqual([AARead('id1', 'actgs')], reads)

    def testFilterRandomSubsetOfZeroFromZeroReads(self):
        """
        It must be possible to select a random subset of zero reads from a set
        of zero reads, where the read count is provided to C{filter} via the
        C{trueLength} argument.
        """
        data = ''
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = FastaReads('filename.fasta')
            result = list(reads.filter(randomSubset=0, trueLength=0))
            self.assertEqual([], result)

    def testFilterRandomSubsetOfTwoFromTwoReads(self):
        """
        It must be possible to select a random subset of two reads from a set
        of two reads, where the read count is provided to C{filter} via the
        C{trueLength} argument.
        """
        data = '\n'.join(['>id1', 'ACGT', '>id2', 'TGCA'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = FastaReads('filename.fasta')
            result = list(reads.filter(randomSubset=2, trueLength=2))
            self.assertEqual([Read('id1', 'ACGT'), Read('id2', 'TGCA')],
                             result)

    def testFilterRandomSubsetOfOneFromTenReads(self):
        """
        It must be possible to select a random subset of one read from a set
        of ten reads, where the read count is provided to C{filter} via the
        C{trueLength} argument.
        """
        data = '\n'.join(['>id', 'ACGT'] * 10)
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = FastaReads('filename.fasta')
            result = list(reads.filter(randomSubset=1, trueLength=10))
            self.assertEqual(1, len(result))


class TestCombineReads(TestCase):
    """
    Tests for the L{dark.fasta.combineReads} function.
    """
    def testNoneNone(self):
        """
        A C{None} FASTA file name and None sequences results in an empty
        FastaReads instance.
        """
        reads = combineReads(None, None)
        self.assertEqual([], list(reads))

    def testNoneEmpty(self):
        """
        A C{None} FASTA file name and an empty sequences list results in an
        empty FastaReads instance.
        """
        reads = list(combineReads(None, []))
        self.assertEqual([], reads)

    def testFileOnly(self):
        """
        If a FASTA file is given but sequences is C{None}, the resulting
        FastaReads must contain the expected read.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(combineReads('filename.fasta', None))
            self.assertEqual([Read('id1', 'ACGT')], reads)

    def testNoUpperCaseFileOnly(self):
        """
        If upperCase is not passed and a FASTA file is given, the resulting
        FastaReads must contain the expected read, in the original case.
        """
        data = '\n'.join(['>id1', 'AcgT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(combineReads('filename.fasta', None))
            self.assertEqual([Read('id1', 'AcgT')], reads)

    def testUpperCaseFileOnly(self):
        """
        When passing upperCase=True and a FASTA file, the resulting
        FastaReads must have the read sequence in uppper case.
        """
        data = '\n'.join(['>id1', 'acgt'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(combineReads('filename.fasta', None, upperCase=True))
            self.assertEqual([Read('id1', 'ACGT')], reads)

    def testSequencesOnly(self):
        """
        A C{None} FASTA file name and a non-empty sequences list results in a
        FastaReads instance with the expected read.
        """
        reads = list(combineReads(None, ['id ACGTSSS'], readClass=AARead))
        self.assertEqual([AARead('id', 'ACGTSSS')], reads)

    def testNoUpperCaseSequencesOnly(self):
        """
        If upperCase is not passed to combineReads the resulting read
        sequences must have their original case.
        """
        reads = list(combineReads(None, ['id aCGt']))
        self.assertEqual([Read('id', 'aCGt')], reads)

    def testUpperCaseSequencesOnly(self):
        """
        Passing upperCase=True to combineReads must result in read sequences
        being upper cased.
        """
        reads = list(combineReads(None, ['id acgt'], upperCase=True))
        self.assertEqual([Read('id', 'ACGT')], reads)

    def testDefaultReadIdPrefix(self):
        """
        A C{None} FASTA file name and a non-empty sequences list with a
        sequence that has no id results in a FastaReads instance with the
        expected read.
        """
        reads = list(combineReads(None, ['ACGT']))
        self.assertEqual([Read('command-line-read-1', 'ACGT')], reads)

    def testCustomReadIdPrefix(self):
        """
        A C{None} FASTA file name and a non-empty sequences list with a
        sequence that has no id, but with a custom read id prefix, results in a
        FastaReads instance with the expected read.
        """
        reads = list(combineReads(None, ['ACGTSSS'], idPrefix='prefix-',
                     readClass=AARead))
        self.assertEqual([AARead('prefix-1', 'ACGTSSS')], reads)

    def testSpecificReadClass(self):
        """
        A specific read class must result in a FastaReads instance with reads
        of that class, both for reads from a FASTA file and from individually
        specified sequences.
        """
        data = '\n'.join(['>id1', 'ACGT'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(combineReads('filename.fasta', ['ACGT'],
                                      readClass=RNARead))
            self.assertTrue(isinstance(reads[0], RNARead))
            self.assertTrue(isinstance(reads[1], RNARead))

import sys
from six.moves import builtins
from six import assertRaisesRegex
from io import BytesIO
import os

from unittest import TestCase, skipUnless
from Bio import SeqIO, bgzf
from contextlib import contextmanager

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import mockOpen

from dark.reads import Read, AARead, DNARead, RNARead, Reads
from dark.fasta import (dedupFasta, dePrefixAndSuffixFasta, fastaSubtract,
                        FastaReads, FastaFaiReads, combineReads, SqliteIndex)
from dark.utils import StringIO

canTestPyfaidx = sys.platform != 'linux'


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

    def testTwoFiles(self):
        """
        It must be possible to read from two FASTA files.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('file1.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n')
                elif self.count == 1:
                    self.test.assertEqual('file2.fasta', filename)
                    self.count += 1
                    return StringIO('>id2\nCAGT\n')
                else:
                    self.test.fail('We are only supposed to be called twice!')

        sideEffect = SideEffect(self)
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = FastaReads(['file1.fasta', 'file2.fasta'])
            self.assertEqual(
                [
                    DNARead('id1', 'ACTG'),
                    DNARead('id2', 'CAGT'),
                ],
                list(reads))


@skipUnless(canTestPyfaidx, 'pyfaidx tests are skipped on Linux')
class TestFastaFaiReads(TestCase):
    """
    Tests for the L{dark.fasta.FastaFaiReads} class.
    """
    def testMissingKey(self):
        """
        If a non-existent sequence id is looked up, a KeyError must be raised.
        """

        pyfaidxIndex = StringIO()

        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return BytesIO(b'>id1\nACTG\n')
                elif self.count == 1:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n')
                elif self.count == 2:
                    self.count += 1
                    return self.manager
                elif self.count == 3:
                    self.count += 1
                    return StringIO(pyfaidxIndex.getvalue())
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        @contextmanager
        def manager():
            yield pyfaidxIndex

        sideEffect = Open(self, manager()).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            reads = FastaFaiReads('filename.fasta')
            error = "^'id2 not in filename\\.fasta\\.'"
            assertRaisesRegex(self, KeyError, error, reads.__getitem__, 'id2')

    def testOneRead(self):
        """
        It must be possible to access a FASTA file with one read like a dict.
        """

        pyfaidxIndex = StringIO()

        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return BytesIO(b'>id1\nACTG\n')
                elif self.count == 1:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n')
                elif self.count == 2:
                    self.count += 1
                    return self.manager
                elif self.count == 3:
                    self.count += 1
                    return StringIO(pyfaidxIndex.getvalue())
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        @contextmanager
        def manager():
            yield pyfaidxIndex

        sideEffect = Open(self, manager()).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            reads = FastaFaiReads('filename.fasta')
            self.assertEqual(DNARead('id1', 'ACTG'), reads['id1'])
            # Check that the fai index was built correctly.
            self.assertEqual(pyfaidxIndex.getvalue(), 'id1\t4\t5\t4\t5\n')

    def testTwoReads(self):
        """
        It must be possible to access a FASTA file with two reads like a dict.
        """

        pyfaidxIndex = StringIO()

        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return BytesIO(b'>id1\nACTG\n>id2\nAACCTTGG\n')
                elif self.count == 1:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAACCTTGG\n')
                elif self.count == 2:
                    self.count += 1
                    return self.manager
                elif self.count == 3:
                    self.count += 1
                    return StringIO(pyfaidxIndex.getvalue())
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        @contextmanager
        def manager():
            yield pyfaidxIndex

        sideEffect = Open(self, manager()).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            reads = FastaFaiReads('filename.fasta')
            self.assertEqual(DNARead('id1', 'ACTG'), reads['id1'])
            self.assertEqual(DNARead('id2', 'AACCTTGG'), reads['id2'])
            # Check that the fai index was built correctly.
            self.assertEqual(pyfaidxIndex.getvalue(),
                             'id1\t4\t5\t4\t5\nid2\t8\t15\t8\t9\n')


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


class TestSqliteIndex(TestCase):
    """
    Tests for the SqliteIndex class.
    """
    def testAddFilename(self):
        """"
        Test the internal _addFilename method.
        """
        index = SqliteIndex(':memory:')
        self.assertEqual(1, index._addFilename('filename1.fasta'))
        self.assertEqual(2, index._addFilename('filename2.fasta'))
        index.close()

    def testAddDuplicateFilename(self):
        """"
        When _addFilename is called twice with the same name, a ValueError
        must be raised.
        """
        index = SqliteIndex(':memory:')
        self.assertEqual(1, index._addFilename('f.fas'))
        error = "^Duplicate file name: 'f.fas'$"
        assertRaisesRegex(self, ValueError, error, index._addFilename, 'f.fas')

    def testGetNonexistentFilename(self):
        """"
        If the internal _getFilename method is called with a file number that
        has not been added, it must return None.
        """
        index = SqliteIndex(':memory:')
        self.assertEqual(None, index._getFilename(1))
        index.close()

    def testGetFilename(self):
        """"
        The internal _getFilename method must return the expected result.
        """
        index = SqliteIndex(':memory:')
        self.assertEqual(1, index._addFilename('filename.fasta'))
        self.assertEqual('filename.fasta', index._getFilename(1))
        index.close()

    def testGetNonexistentFileNumber(self):
        """"
        If the internal _getFileNumber method is called with a file whose name
        has not been added, it must return None.
        """
        index = SqliteIndex(':memory:')
        self.assertEqual(None, index._getFileNumber('filename.fasta'))
        index.close()

    def testGetFileNumber(self):
        """"
        The internal _getFileNumber method must return the expected result.
        """
        index = SqliteIndex(':memory:')
        self.assertEqual(1, index._addFilename('filename.fasta'))
        self.assertEqual(1, index._getFileNumber('filename.fasta'))
        index.close()

    def testBZ2File(self):
        """"
        Trying to add a .bz2 file must result in a ValueError.
        """
        index = SqliteIndex(':memory:')
        error = ('^Compressed FASTA is only supported in BGZF format\\. Use '
                 'bgzip to compresss your FASTA\\.$')
        assertRaisesRegex(self, ValueError, error, index.addFile, 'file.bz2')

    def testAddOneFile(self):
        """"
        Test the creation of an index with sequences added from one file.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAACCTTGG\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            self.assertEqual(2, index.addFile('filename.fasta'))
            index.close()

    def testAddFileWithDuplicateSequence(self):
        """"
        If a sequence id is duplicated in a FASTA file, a ValueError must be
        raised.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id1\nAACCTTGG\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            error = ("^FASTA sequence id 'id1' found twice in file "
                     "'filename.fasta'\\.$")
            assertRaisesRegex(self, ValueError, error, index.addFile,
                              'filename.fasta')
            index.close()

    def testAddFilesWithDuplicateSequence(self):
        """"
        If a sequence id occurs in more than one FASTA file, a ValueError must
        be raised.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename1.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAACCTTGG\n')
                elif self.count == 1:
                    self.test.assertEqual('filename2.fasta', filename)
                    self.count += 1
                    return StringIO('>id2\nAAACCC\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename1.fasta')
            error = ("^FASTA sequence id 'id2', found in file "
                     "'filename2\\.fasta', was previously added from file "
                     "'filename1\\.fasta'\\.$")
            assertRaisesRegex(self, ValueError, error, index.addFile,
                              'filename2.fasta')
            index.close()

    def testAddDuplicateFile(self):
        """"
        If a filename is passed to addFile more than once, a ValueError must
        be raised.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAACCTTGG\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            self.assertEqual(2, index.addFile('filename.fasta'))
            error = "^Duplicate file name: 'filename\\.fasta'$"
            assertRaisesRegex(self, ValueError, error, index._addFilename,
                              'filename.fasta')
            index.close()

    def testFind(self):
        """"
        The _find method must return the expected filename and offset.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAACCTTGG\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename.fasta')
            self.assertEqual(('filename.fasta', 5), index._find('id1'))
            self.assertEqual(('filename.fasta', 15), index._find('id2'))
            index.close()

    def testFindWithTwoFiles(self):
        """"
        The _find method must return the expected filename and offset when
        sequences are added from two files.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename1.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAACCTTGG\n')
                elif self.count == 1:
                    self.test.assertEqual('filename2.fasta', filename)
                    self.count += 1
                    return StringIO('>sequence3\nAAACCC\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename1.fasta')
            index.addFile('filename2.fasta')
            self.assertEqual(('filename1.fasta', 5), index._find('id1'))
            self.assertEqual(('filename1.fasta', 15), index._find('id2'))
            self.assertEqual(('filename2.fasta', 11), index._find('sequence3'))
            index.close()

    def testDictLookupSequenceCrossesNewlines(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected read when the sequence spans multiple lines of the input file,
        including lines ending in \n and \r\n.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0 or self.count == 1:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\r\nCCCC\nGGG\n>id2\nAACCTG\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename.fasta')
            self.assertEqual(DNARead('id1', 'ACTGCCCCGGG'), index['id1'])
            index.close()

    def testDictLookupWithFastaDirectory(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected read, obtained from the expected file name, when a FASTA base
        directory is specified.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('/tmp/f.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\r\nCCCC\nGGG\n>id2\nAACCTG\n')
                if self.count == 1:
                    self.test.assertEqual(
                        os.path.join('/usr/local/fasta', 'f.fasta'), filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\r\nCCCC\nGGG\n>id2\nAACCTG\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:', fastaDirectory='/usr/local/fasta')
            index.addFile('/tmp/f.fasta')
            self.assertEqual(DNARead('id1', 'ACTGCCCCGGG'), index['id1'])
            index.close()

    def testDictLookupSequenceLastInFile(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected read when the sequence spans multiple lines and is the last
        one in the input file.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0 or self.count == 1:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\r\nCCCC\n>id2\nAACCTG\nAAA\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename.fasta')
            self.assertEqual(DNARead('id2', 'AACCTGAAA'), index['id2'])
            index.close()

    def testDictLookupSequenceMiddleOfThree(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected read when the sequence spans multiple lines and is the middle
        one of three sequences in the input file.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0 or self.count == 1:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO(
                        '>id1\nACTG\nCCCC\n>id2\nAACCTG\nAAA\n>id3\nAAA\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename.fasta')
            self.assertEqual(DNARead('id2', 'AACCTGAAA'), index['id2'])
            index.close()

    def testDictLookupWithTwoFiles(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected reads when sequences are added from two files.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0 or self.count == 2 or self.count == 3:
                    self.test.assertEqual('filename1.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAACCTTGG\n')
                elif self.count == 1 or self.count == 4:
                    self.test.assertEqual('filename2.fasta', filename)
                    self.count += 1
                    return StringIO('>seq3\nAAACCC\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename1.fasta')
            index.addFile('filename2.fasta')
            self.assertEqual(DNARead('id1', 'ACTG'), index['id1'])
            self.assertEqual(DNARead('id2', 'AACCTTGG'), index['id2'])
            self.assertEqual(DNARead('seq3', 'AAACCC'), index['seq3'])
            index.close()

    def testDictLookupSpecificReadClass(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected read type.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0 or self.count == 1:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nMM\n>id2\n')
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:', readClass=AARead)
            index.addFile('filename.fasta')
            result = index['id1']
            self.assertTrue(isinstance(result, AARead))
            self.assertEqual(AARead('id1', 'MM'), result)
            index.close()

    def testDictLookupGzipDataWithBGZsuffix(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected read when the index file is in BGZF format and has a .bgz
        suffix.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count <= 1:
                    self.test.assertEqual('filename.fasta.bgz', filename)
                    self.count += 1
                    writerIO = BytesIO()
                    writer = bgzf.BgzfWriter(fileobj=writerIO)
                    writer.write(b'>id0\nAC\n')
                    writer.flush()
                    fileobj = BytesIO(writerIO.getvalue())
                    fileobj.mode = 'rb'
                    return bgzf.BgzfReader(fileobj=fileobj)
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(bgzf, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename.fasta.bgz')
            self.assertEqual(DNARead('id0', 'AC'), index['id0'])
            index.close()

    def testDictLookupGzipData(self):
        """"
        The __getitem__ method (i.e., dictionary-like lookup) must return the
        expected reads when sequences span multiple lines of the input file,
        and include lines ending in \n and \r\n and have been compressed with
        bgzip, including when sequences are more than 64K bytes into the input
        file.
        """
        class Open(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count <= 4:
                    self.test.assertEqual('filename.fasta.gz', filename)
                    self.count += 1
                    writerIO = BytesIO()
                    writer = bgzf.BgzfWriter(fileobj=writerIO)
                    writer.write(
                        b'>id0\nAC\n' +
                        b'>id1\n' + (b'A' * 70000) + b'\n' +
                        b'>id2\r\nACTG\r\nCCCC\r\nGGG\r\n' +
                        b'>id3\nAACCTG\n')
                    writer.flush()
                    fileobj = BytesIO(writerIO.getvalue())
                    fileobj.mode = 'rb'
                    return bgzf.BgzfReader(fileobj=fileobj)
                else:
                    self.test.fail(
                        'Open called too many times. Filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        sideEffect = Open(self).sideEffect
        with patch.object(bgzf, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            index = SqliteIndex(':memory:')
            index.addFile('filename.fasta.gz')
            self.assertEqual(DNARead('id0', 'AC'), index['id0'])
            self.assertEqual(DNARead('id1', 'A' * 70000), index['id1'])
            self.assertEqual(DNARead('id2', 'ACTGCCCCGGG'), index['id2'])
            self.assertEqual(DNARead('id3', 'AACCTG'), index['id3'])
            index.close()

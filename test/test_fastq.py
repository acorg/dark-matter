import bz2

from six.moves import builtins

from dark.reads import AARead, DNARead, RNARead
from dark.fastq import FastqReads

from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import mockOpen


class BZ2(object):
    """
    A BZ2File mock.
    """
    def __init__(self, data):
        self._data = data
        self._index = 0

    def readline(self):
        self._index += 1
        try:
            return self._data[self._index - 1]
        except IndexError:
            return None

    def __iter__(self):
        index = self._index
        self._index = len(self._data)
        return iter(self._data[index:])


class TestFastqReads(TestCase):
    """
    Tests for the L{dark.fastq.FastqReads} class.
    """

    def testEmpty(self):
        """
        An empty FASTQ file results in an empty iterator.
        """
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            reads = FastqReads('filename.fastq')
            self.assertEqual([], list(reads))

    def testOneRead(self):
        """
        A FASTQ file with one read must be read properly.
        """
        data = '\n'.join(['@id1', 'ACGT', '+', '!!!!'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastqReads('filename.fastq'))
            self.assertEqual([DNARead('id1', 'ACGT', '!!!!')], reads)

    def testTwoReads(self):
        """
        A FASTQ file with two reads must be read properly and its
        sequences must be returned in the correct order.
        """
        data = '\n'.join(['@id1', 'ACGT', '+', '!!!!',
                          '@id2', 'TGCA', '+', '????'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastqReads('filename.fastq'))
            self.assertEqual(2, len(reads))
            self.assertEqual([DNARead('id1', 'ACGT', '!!!!'),
                              DNARead('id2', 'TGCA', '????')], reads)

    def testTypeDefaultsToDNA(self):
        """
        A FASTQ file whose type is not specified must result in reads that
        are instances of DNARead.
        """
        data = '\n'.join(['@id1', 'ACGT', '+', '!!!!'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastqReads('filename.fastq'))
            self.assertTrue(isinstance(reads[0], DNARead))

    def testTypeAA(self):
        """
        A FASTQ file whose read class is AARead must result in reads that
        are instances of AARead.
        """
        data = '\n'.join(['@id1', 'ACGT', '+', '!!!!'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastqReads('filename.fastq', AARead))
            self.assertTrue(isinstance(reads[0], AARead))

    def testTypeDNA(self):
        """
        A FASTQ file whose read class is DNARead must result in reads that
        are instances of DNARead.
        """
        data = '\n'.join(['@id1', 'ACGT', '+', '!!!!'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastqReads('filename.fastq', DNARead))
            self.assertTrue(isinstance(reads[0], DNARead))

    def testTypeRNA(self):
        """
        A FASTQ file whose read class is RNARead must result in reads that
        are instances of RNARead.
        """
        data = '\n'.join(['@id1', 'ACGT', '+', '!!!!'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = list(FastqReads('filename.fastq', RNARead))
            self.assertTrue(isinstance(reads[0], RNARead))

    def testTwoFiles(self):
        """
        It must be possible to read from two FASTQ files.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename):
                if self.count == 0:
                    self.test.assertEqual('file1.fastq.bz2', filename)
                    self.count += 1
                    return BZ2(['@id1\n', 'ACTG\n', '+\n', '!!!!\n'])
                elif self.count == 1:
                    self.test.assertEqual('file2.fastq.bz2', filename)
                    self.count += 1
                    return BZ2(['@id2\n', 'CAGT\n', '+\n', '!!!!\n'])
                else:
                    self.fail('We are only supposed to be called twice!')

        sideEffect = SideEffect(self)
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = FastqReads(['file1.fastq.bz2', 'file2.fastq.bz2'])
            self.assertEqual(
                [
                    DNARead('id1', 'ACTG', '!!!!'),
                    DNARead('id2', 'CAGT', '!!!!'),
                ],
                list(reads))

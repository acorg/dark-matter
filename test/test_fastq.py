import builtins
from io import BytesIO
from unittest import TestCase
from unittest.mock import mock_open, patch

from dark.fastq import FastqReads
from dark.reads import AARead, DNARead, RNARead


class TestFastqReads(TestCase):
    """
    Tests for the L{dark.fastq.FastqReads} class.
    """

    def testEmpty(self):
        """
        An empty FASTQ file results in an empty iterator.
        """
        with patch.object(builtins, "open", mock_open(read_data=b"")):
            reads = FastqReads("filename.fastq")
            self.assertEqual([], list(reads))

    def testOneRead(self):
        """
        A FASTQ file with one read must be read properly.
        """
        data = b"\n".join([b"@id1", b"ACGT", b"+", b"!!!!"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(FastqReads("filename.fastq"))
            self.assertEqual([DNARead("id1", "ACGT", "!!!!")], reads)

    def testTwoReads(self):
        """
        A FASTQ file with two reads must be read properly and its
        sequences must be returned in the correct order.
        """
        data = b"\n".join(
            [b"@id1", b"ACGT", b"+", b"!!!!", b"@id2", b"TGCA", b"+", b"????"]
        )
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(FastqReads("filename.fastq"))
            self.assertEqual(2, len(reads))
            self.assertEqual(
                [DNARead("id1", "ACGT", "!!!!"), DNARead("id2", "TGCA", "????")], reads
            )

    def testTypeDefaultsToDNA(self):
        """
        A FASTQ file whose type is not specified must result in reads that
        are instances of DNARead.
        """
        data = b"\n".join([b"@id1", b"ACGT", b"+", b"!!!!"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(FastqReads("filename.fastq"))
            self.assertTrue(isinstance(reads[0], DNARead))

    def testTypeAA(self):
        """
        A FASTQ file whose read class is AARead must result in reads that
        are instances of AARead.
        """
        data = b"\n".join([b"@id1", b"ACGT", b"+", b"!!!!"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(FastqReads("filename.fastq", AARead))
            self.assertTrue(isinstance(reads[0], AARead))

    def testTypeDNA(self):
        """
        A FASTQ file whose read class is DNARead must result in reads that
        are instances of DNARead.
        """
        data = b"\n".join([b"@id1", b"ACGT", b"+", b"!!!!"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(FastqReads("filename.fastq", DNARead))
            self.assertTrue(isinstance(reads[0], DNARead))

    def testTypeRNA(self):
        """
        A FASTQ file whose read class is RNARead must result in reads that
        are instances of RNARead.
        """
        data = b"\n".join([b"@id1", b"ACGT", b"+", b"!!!!"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(FastqReads("filename.fastq", RNARead))
            self.assertTrue(isinstance(reads[0], RNARead))

    def testTwoFiles(self):
        """
        It must be possible to read from two FASTQ files.
        """

        class SideEffect:
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, mode):
                assert "b" in mode
                if self.count == 0:
                    self.test.assertEqual("file1.fastq", filename)
                    self.count += 1
                    return BytesIO(b"@id1\nACTG\n+\n!!!!\n")
                elif self.count == 1:
                    self.test.assertEqual("file2.fastq", filename)
                    self.count += 1
                    return BytesIO(b"@id2\nCAGT\n+\n!!!!\n")
                else:
                    self.test.fail("We are only supposed to be called twice!")

        sideEffect = SideEffect(self)
        with patch.object(builtins, "open") as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = FastqReads(["file1.fastq", "file2.fastq"])
            self.assertEqual(
                [
                    DNARead("id1", "ACTG", "!!!!"),
                    DNARead("id2", "CAGT", "!!!!"),
                ],
                list(reads),
            )

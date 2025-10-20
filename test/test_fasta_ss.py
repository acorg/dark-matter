import builtins
from io import BytesIO
from unittest import TestCase
from unittest.mock import mock_open, patch

from dark.fasta_ss import SSFastaReads
from dark.reads import SSAARead


class TestSSFastaReads(TestCase):
    """
    Tests for the L{dark.fasta.SSFastaReads} class.
    """

    def testEmpty(self):
        """
        An empty PDB FASTA file results in an empty iterator.
        """
        data = b""
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = SSFastaReads("x.fasta")
            self.assertEqual([], list(reads))

    def testOddNumberOfRecords(self):
        """
        Trying to parse a PDB FASTA file with an odd number of records must
        raise a ValueError.
        """
        data = b"\n".join([b">seq1", b"REDD", b">str1", b"HH--", b">seq2", b"REAA"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            error = "^Structure file 'x.fasta' has an odd number of records\\.$"
            self.assertRaisesRegex(ValueError, error, list, SSFastaReads("x.fasta"))

    def testUnequalSequenceAndStructureLengths(self):
        """
        Trying to parse a PDB FASTA file that has a sequence whose structure
        is of a different length must raise a ValueError.
        """
        data = b"\n".join(
            [b">seq1", b"REDD", b">str1", b"HH--", b">seq2", b"REAA", b">str2", b"HH"]
        )
        with patch.object(builtins, "open", mock_open(read_data=data)):
            error = (
                "Sequence 'seq2' length \\(4\\) is not equal to structure "
                "'str2' length \\(2\\) in input file 'x\\.fasta'\\.$"
            )
            self.assertRaisesRegex(ValueError, error, list, SSFastaReads("x.fasta"))

    def testOneRead(self):
        """
        A PDB FASTA file with one read must be read properly.
        """
        data = b"\n".join([b">seq1", b"REDD", b">str1", b"HH--"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(SSFastaReads("x.fasta"))
            self.assertEqual([SSAARead("seq1", "REDD", "HH--")], reads)

    def testNoQuality(self):
        """
        A PDB FASTA file read must not have any quality information.
        """
        data = b"\n".join([b">seq1", b"REDD", b">str1", b"HH--"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(SSFastaReads("x.fasta"))
            self.assertIs(None, reads[0].quality)

    def testTwoReads(self):
        """
        A PDB FASTA file with two reads must be read properly and its
        sequences must be returned in the correct order.
        """
        data = b"\n".join(
            [b">seq1", b"REDD", b">str1", b"HH--", b">seq2", b"REAA", b">str2", b"HHEE"]
        )
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(SSFastaReads("x.fasta"))
            self.assertEqual(2, len(reads))
            self.assertEqual(
                [SSAARead("seq1", "REDD", "HH--"), SSAARead("seq2", "REAA", "HHEE")],
                reads,
            )

    def testTypeDefaultsToSSAARead(self):
        """
        A PDB FASTA file whose type is not specified must result in reads that
        are instances of SSAARead.
        """
        data = b"\n".join([b">seq1", b"REDD", b">str1", b"HH--"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(SSFastaReads("x.fasta"))
            self.assertTrue(isinstance(reads[0], SSAARead))

    def testReadClass(self):
        """
        A PDB FASTA file whose read class is something other than SSAARead must
        result in reads that are instances of that class.
        """

        class ReadClass:
            def __init__(self, id, sequence, structure):
                pass

        data = b"\n".join([b">seq1", b"RRRR", b">str1", b"HHHH"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(SSFastaReads("x.fasta", readClass=ReadClass))
            self.assertTrue(isinstance(reads[0], ReadClass))

    def testConvertLowerToUpperCaseIfSpecified(self):
        """
        A read sequence and structure must be converted from lower to upper
        case if requested.
        """
        data = b"\n".join([b">seq1", b"rrrff", b">str1", b"hheeh"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(SSFastaReads("x.fasta", upperCase=True))
            self.assertEqual([SSAARead("seq1", "RRRFF", "HHEEH")], reads)

    def testDontConvertLowerToUpperCaseIfNotSpecified(self):
        """
        A read sequence and its structure must not be converted from lower to
        upper case if the conversion is not requested.
        """
        data = b"\n".join([b">seq1", b"rrFF", b">str1", b"HHee"])
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reads = list(SSFastaReads("x.fasta"))
            self.assertEqual([SSAARead("seq1", "rrFF", "HHee")], reads)

    def testTwoFiles(self):
        """
        It must be possible to read from two FASTA files.
        """

        class SideEffect:
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, **kwargs):
                if self.count == 0:
                    self.test.assertEqual("file1.fasta", filename)
                    self.count += 1
                    return BytesIO(b">id1\nACTG\n>id1\nhhhh\n")
                elif self.count == 1:
                    self.test.assertEqual("file2.fasta", filename)
                    self.count += 1
                    return BytesIO(b">id2\nCAGT\n>id2\neeee\n")
                else:
                    self.test.fail("We are only supposed to be called twice!")

        sideEffect = SideEffect(self)
        with patch.object(builtins, "open") as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = SSFastaReads(["file1.fasta", "file2.fasta"])
            self.assertEqual(
                [
                    SSAARead("id1", "ACTG", "hhhh"),
                    SSAARead("id2", "CAGT", "eeee"),
                ],
                list(reads),
            )

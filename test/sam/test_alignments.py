from unittest import TestCase
from mock import patch
from ..mocking import mockOpen
import sample_data
import sample_data_params
import sample_data_1

from dark.reads import Read, Reads
from dark.sam.alignments import SamReadsAlignments


class TestSamReadsAlignments(TestCase):
    """
    Test the SamReadsAlignments class.
    """

    def testEmptySAMInput(self):
        """
        When a SAM input file is empty, a C{ValueError} must be raised
        on trying to read it.
        """
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            error = "SAM file 'file.SAM' was empty."
            self.assertRaisesRegexp(
                ValueError, error, SamReadsAlignments, reads, 'file.SAM')

    def testNonSAMInput(self):
        """
        When given a file whose contents are not in SAM file format,
        attempting to read the BLAST hits from it must raise a
        C{ValueError}.
        """
        mockOpener = mockOpen(read_data='not SAM\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            error = ("No header lines in file.SAM")
            self.assertRaisesRegexp(
                ValueError, error, SamReadsAlignments, reads, 'file.SAM')

    def testApplicationParams(self):
        """
        Alignment parameters must be extracted from the input SAM file and stored
        correctly.
        """
        mockOpener = mockOpen(read_data=open(sample_data))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = SamReadsAlignments(reads, 'file.SAM')
            params = {'VN': '1.3', 'SO': 'unsorted', 'GO': 'none',
                      'gi|887699|gb|DQ37780': 37000,
                      'gi|639163157|ref|NC_024124.1|': 35000, 'app': 'bfast',
                      'VN': '0.7.0a'}
            self.assertEqual(params, readsAlignments.params.applicationParams)

    def testSAMParamsButNoHits(self):
        """
        When alignment parameters are present in the input but there are no
        records, the __iter__ method of a L{SamReadsAlignments} instance must
        not yield anything.
        """
        mockOpener = mockOpen(read_data=open(sample_data_params))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = SamReadsAlignments(reads, 'file.SAM')
            self.assertEqual([], list(readsAlignments))

    def testTooManyReads(self):
        """
        If a SAM file contains a parameters section and one hit, but there
        is more than one read, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(read_data=open(sample_data_1))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'G' * 70))
            error = ("Reads iterator contained more reads than the number of "
                     "BLAST records found \(1\)\. First unknown read id is "
                     "'id1'\.")
            readsAlignments = SamReadsAlignments(reads, 'file.SAM')
            self.assertRaisesRegexp(ValueError, error, list, readsAlignments)

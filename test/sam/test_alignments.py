from unittest import TestCase
from mock import patch
from ..mocking import mockOpen

from dark.reads import Read, Reads
from dark.sam.alignments import SamReadsAlignments

SAMPLE_DATA = """\
@HD VN:1.3  SO:unsorted GO:none
@SQ SN:gi|887699|gb|DQ37780 LN:37000
@SQ SN:gi|639163157|ref|NC_024124.1| LN:35000
@PG ID:bfast    VN:0.7.0a
id0 0 gi|887699|gb|DQ37780 15363 5 5M2I9M1I7M1D6M * 0 \
TACCCTGCGGCCCGCTACGGCTGGTCTCCA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MD:Z:5^TG8T^A4T8
id1 0 gi|639163157|ref|NC_024124.1| 11362 13 5M2I9M1I7M1D6M * 0 \
TACCCTGCGGCCCGCTACGGCTGGTCTCCA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MD:Z:5^TG8T^A4T8
id1 0 gi|887699|gb|DQ37780 11362 13 5M2I9M1I7M1D6M * 0 \
TACCCTGCGGCCCGCTACGGCTGGTCTCCA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MD:Z:5^TG8T^A4T8
"""

SAMPLE_DATA_1 = """\
@HD VN:1.3  SO:unsorted GO:none
@SQ SN:gi|887699|gb|DQ37780 LN:37000
@PG ID:bfast    VN:0.7.0a
id0 0 gi|887699|gb|DQ37780 15363 5 5M2I9M1I7M1D6M * 0 \
TACCCTGCGGCCCGCTACGGCTGGTCTCCA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MD:Z:5^TG8T^A4T8
"""

SAMPLE_DATA_PARAMS = """\
@HD VN:1.3  SO:unsorted GO:none
@SQ SN:gi|887699|gb|DQ37780 LN:37000
@SQ SN:gi|639163157|ref|NC_024124.1| LN:35000
@PG ID:bfast    VN:0.7.0a
"""


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
            error = "SAM file file.SAM was empty."
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
            error = "SAM file file.SAM does not contain at least 11 fields."
            self.assertRaisesRegexp(
                AssertionError, error, SamReadsAlignments, reads, 'file.SAM')

    def testApplicationParams(self):
        """
        Alignment parameters must be extracted from the input SAM file and
        stored correctly.
        """

        mockOpener = mockOpen(read_data=SAMPLE_DATA)
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = SamReadsAlignments(reads, 'file.SAM')
            params = {'VN': '1.3', 'SO': 'unsorted', 'GO': 'none',
                      'gi|887699|gb|DQ37780': 37000,
                      'gi|639163157|ref|NC_024124.1|': 35000,
                      'application': 'bfast', 'VN': '0.7.0a'}
            self.assertEqual(params, readsAlignments.params.applicationParams)

    def testSAMParamsButNoHits(self):
        """
        When alignment parameters are present in the input but there are no
        records, the __iter__ method of a L{SamReadsAlignments} instance must
        not yield anything.
        """
        mockOpener = mockOpen(read_data=SAMPLE_DATA_PARAMS)
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = SamReadsAlignments(reads, 'file.SAM')
            self.assertEqual([], list(readsAlignments))

    def testTooManyReads(self):
        """
        If a SAM file contains a parameters section and one hit, but there
        is more than one read, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(read_data=SAMPLE_DATA_1)
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'G' * 70))
            error = ("Reads iterator contained more reads than the number of "
                     "BLAST records found \(1\)\. First unknown read id is "
                     "'id1'\.")
            readsAlignments = SamReadsAlignments(reads, 'file.SAM')
            self.assertRaisesRegexp(ValueError, error, list, readsAlignments)

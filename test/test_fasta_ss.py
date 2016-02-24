import six
from six.moves import builtins
from six import StringIO
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import mockOpen

from dark.reads import SSAARead
from dark.fasta_ss import SSFastaReads


class TestSSFastaReads(TestCase):
    """
    Tests for the L{dark.fasta.SSFastaReads} class.
    """

    def testEmpty(self):
        """
        An empty PDB FASTA file results in an empty iterator.
        """
        data = StringIO()
        reads = SSFastaReads(data)
        self.assertEqual([], list(reads))

    def testOddNumberOfRecords(self):
        """
        Trying to parse a PDB FASTA file with an odd number of records must
        raise a ValueError.
        """
        data = '\n'.join(['>seq1', 'REDD', '>str1', 'HH--',
                          '>seq2', 'REAA'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = "^Structure file 'x.fasta' has an odd number of records\.$"
            six.assertRaisesRegex(self, ValueError, error, list,
                                  SSFastaReads('x.fasta'))

    def testUnequalSequenceAndStructureLengths(self):
        """
        Trying to parse a PDB FASTA file that has a sequence whose structure
        is of a different length must raise a ValueError.
        """
        data = '\n'.join(['>seq1', 'REDD', '>str1', 'HH--',
                          '>seq2', 'REAA', '>str2', 'HH'])
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            error = ("Sequence 'seq2' length \(4\) is not equal to structure "
                     "'str2' length \(2\) in input file 'x\.fasta'\.$")
            six.assertRaisesRegex(self, ValueError, error, list,
                                  SSFastaReads('x.fasta'))

    def testOneRead(self):
        """
        A PDB FASTA file with one read must be read properly.
        """
        data = StringIO('\n'.join(['>seq1', 'REDD', '>str1', 'HH--']))
        reads = list(SSFastaReads(data))
        self.assertEqual([SSAARead('seq1', 'REDD', 'HH--')], reads)

    def testNoQuality(self):
        """
        A PDB FASTA file read must not have any quality information.
        """
        data = StringIO('\n'.join(['>seq1', 'REDD', '>str1', 'HH--']))
        reads = list(SSFastaReads(data))
        self.assertIs(None, reads[0].quality)

    def testTwoReads(self):
        """
        A PDB FASTA file with two reads must be read properly and its
        sequences must be returned in the correct order.
        """
        data = StringIO('\n'.join(['>seq1', 'REDD', '>str1', 'HH--',
                                   '>seq2', 'REAA', '>str2', 'HHEE']))
        reads = list(SSFastaReads(data))
        self.assertEqual(2, len(reads))
        self.assertEqual([SSAARead('seq1', 'REDD', 'HH--'),
                          SSAARead('seq2', 'REAA', 'HHEE')],
                         reads)

    def testTypeDefaultsToSSAARead(self):
        """
        A PDB FASTA file whose type is not specified must result in reads that
        are instances of SSAARead.
        """
        data = StringIO('\n'.join(['>seq1', 'REDD', '>str1', 'HH--']))
        reads = list(SSFastaReads(data))
        self.assertTrue(isinstance(reads[0], SSAARead))

    def testReadClass(self):
        """
        A PDB FASTA file whose read class is something other than SSAARead must
        result in reads that are instances of that class.
        """
        class ReadClass:
            def __init__(self, id, sequence, structure):
                pass

        data = StringIO('\n'.join(['>seq1', 'RRRR', '>str1', 'HHHH']))
        reads = list(SSFastaReads(data, readClass=ReadClass, checkAlphabet=0))
        self.assertTrue(isinstance(reads[0], ReadClass))

    def testAlphabetIsCheckedAndRaisesValueErrorOnFirstRead(self):
        """
        The default behavior of a SSFastaReads instance is to check to ensure
        its sequences have the correct alphabet and to raise ValueError if not.
        A non-alphabetic character in the first read must be detected.
        """
        data = StringIO('\n'.join(['>seq1', 'at-at', '>str1', 'HH-HH']))
        error = ("^Read alphabet \('-AT'\) is not a subset of expected "
                 "alphabet \('ACDEFGHIKLMNPQRSTVWY'\) for read class "
                 "SSAARead\.$")
        six.assertRaisesRegex(self, ValueError, error, list,
                              SSFastaReads(data))

    def testAlphabetIsCheckedAndRaisesValueErrorOnSecondRead(self):
        """
        The default behavior of a SSFastaReads instance is to check to ensure
        its sequences have the correct alphabet and to raise ValueError if not.
        A non-alphabetic character in the second read must be detected.
        """
        data = StringIO('\n'.join(['>seq1', 'rrrr', '>str1', 'hhhh',
                                   '>seq2', 'a-at', '>str2', 'hhhh']))
        error = ("^Read alphabet \('-AT'\) is not a subset of expected "
                 "alphabet \('ACDEFGHIKLMNPQRSTVWY'\) for read class "
                 "SSAARead\.$")
        six.assertRaisesRegex(self, ValueError, error, list,
                              SSFastaReads(data))

    def testDisableAlphabetChecking(self):
        """
        It must be possible to have a SSFastaReads instance not do alphabet
        checking, if requested (by passing checkAlphabet=0).
        """
        data = StringIO('\n'.join(['>seq1', 'rr-rr', '>str1', 'hh-hh']))
        self.assertEqual(1, len(list(SSFastaReads(data, checkAlphabet=0))))

    def testOnlyCheckSomeAlphabets(self):
        """
        It must be possible to have the alphabets of only a certain number of
        reads checked. A non-alphabetic character in a later read must not
        stop that read from being processed.
        """
        data = StringIO('\n'.join(['>seq1', 'rrrr', '>str1', 'hhhh',
                                   '>seq2', 'r-rr', '>str2', 'h-hh']))
        reads = list(SSFastaReads(data, checkAlphabet=1))
        self.assertEqual(2, len(reads))
        self.assertEqual('r-rr', reads[1].sequence)

    def testConvertLowerToUpperCaseIfSpecified(self):
        """
        A read sequence and structure must be converted from lower to upper
        case if requested.
        """
        data = StringIO('\n'.join(['>seq1', 'rrrff', '>str1', 'hheeh']))
        reads = list(SSFastaReads(data, upperCase=True))
        self.assertEqual([SSAARead('seq1', 'RRRFF', 'HHEEH')], reads)

    def testDontConvertLowerToUpperCaseIfNotSpecified(self):
        """
        A read sequence and its structure must not be converted from lower to
        upper case if the conversion is not requested.
        """
        data = StringIO('\n'.join(['>seq1', 'rrFF', '>str1', 'HHee']))
        reads = list(SSFastaReads(data))
        self.assertEqual([SSAARead('seq1', 'rrFF', 'HHee')], reads)

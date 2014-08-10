from unittest import TestCase
from mock import patch, call

from mocking import mockOpen
from dark.reads import Read, Reads


class TestRead(TestCase):
    """
    Test the Read class.
    """

    def testUnknownType(self):
        """
        Attempting to construct a read with an unknown type must raise a
        ValueError.
        """
        error = "Unknown sequence type 'weird'"
        with self.assertRaisesRegexp(ValueError, error):
            Read('id', 'ACGT', '!!!!', type='weird')

    def testUnequalLengths(self):
        """
        Attempting to construct a read whose sequence and quality strings are
        of different lengths must raise a ValueError.
        """
        error = 'Invalid read: sequence length \(4\) != quality length \(3\)'
        with self.assertRaisesRegexp(ValueError, error):
            Read('id', 'ACGT', '!!!')

    def testNoQuality(self):
        """
        If no quality information is given, the read's 'quality' attribute must
        be None.
        """
        read = Read('id', 'ACGT')
        self.assertIs(None, read.quality)

    def testNoType(self):
        """
        If no type is given, the read's type must be DNA.
        """
        read = Read('id', 'ACGT')
        self.assertIs('dna', read.type)

    def testConvertToUpperCase(self):
        """
        The sequence passed to Read must be converted to upper case.
        """
        read = Read('id', 'acgt')
        self.assertEqual('ACGT', read.sequence)

    def testExpectedAttributes(self):
        """
        After constructing a read, the expected attributes must be present.
        """
        read = Read('id', 'ACGT', '!!!!')
        self.assertEqual('id', read.id)
        self.assertEqual('ACGT', read.sequence)
        self.assertEqual('!!!!', read.quality)
        self.assertEqual('dna', read.type)

    def testLength(self):
        """
        len() must return the length of a read's sequence.
        """
        read = Read('id', 'ACGT', '!!!!')
        self.assertEqual(4, len(read))

    def testEqualityWithDifferingIds(self):
        """
        If two Read instances have different ids, they should not be
        considered equal.
        """
        self.assertNotEqual(Read('id1', 'AC'), Read('id2', 'AC'))

    def testEqualityWithDifferingSequences(self):
        """
        If two Read instances have different sequences, they should not be
        considered equal.
        """
        self.assertNotEqual(Read('id1', 'AA'), Read('id1', 'CC'))

    def testEqualityWithDifferingQuality(self):
        """
        If two Read instances have different quality, they should not be
        considered equal.
        """
        self.assertNotEqual(Read('id1', 'AC', 'qq'), Read('id1', 'AC', 'rr'))

    def testEqualityWithOneOmittedQuality(self):
        """
        If two Read instances have different quality (one is omitted), they
        should not be considered equal.
        """
        self.assertNotEqual(Read('id1', 'AC'), Read('id1', 'AC', 'rr'))

    def testEqualityWithNoQuality(self):
        """
        If two Read instances have the same id and sequence and neither has a
        quality, they should be considered equal.
        """
        self.assertEqual(Read('id1', 'AC'), Read('id1', 'AC'))

    def testEquality(self):
        """
        If two Read instances have the same id, sequence, and quality, they
        should be considered equal.
        """
        self.assertEqual(Read('id1', 'AC', 'qq'), Read('id1', 'AC', 'qq'))

    def testReverseComplementAA(self):
        """
        The reverseComplement function must raise a C{ValueError} when called
        on an amino acid sequence.
        """
        read = Read('id', 'atcg', type='aa')
        error = 'Cannot reverse complement an amino acid sequence'
        with self.assertRaisesRegexp(ValueError, error):
            read.reverseComplement()

    def testReverseComplementReversesQuality(self):
        """
        The reverseComplement function must return a reversed quality string.
        """
        read = Read('id', 'atcg', quality='!@#$')
        self.assertEqual('$#@!', read.reverseComplement().quality)

    def testReverseComplementDNA(self):
        """
        The reverseComplement function must work for DNA
        """
        read = Read('id', 'atcg', quality='!@#$', type='dna')
        self.assertEqual('CGAT', read.reverseComplement().sequence)

    def testReverseComplementAmbiguousDNA(self):
        """
        The reverseComplement function must work for DNA that includes
        ambiguous bases.
        """
        read = Read('id', 'atcgmrwsvhxn', type='dna')
        self.assertEqual('NXDBSWYKCGAT', read.reverseComplement().sequence)

    def testReverseComplementRNA(self):
        """
        The reverseComplement function must work for RNA
        """
        read = Read('id', 'aucg', type='rna')
        self.assertEqual('CGAU', read.reverseComplement().sequence)

    def testReverseComplementAmbiguousRNA(self):
        """
        The reverseComplement function must work for RNA that includes
        ambiguous bases.
        """
        read = Read('id', 'aucgmrwsykvhxn', type='rna')
        self.assertEqual('NXDBMRSWYKCGAU', read.reverseComplement().sequence)


class TestReads(TestCase):
    """
    Test the Reads class.
    """

    def testNoReads(self):
        """
        A Reads instance with no reads should return an empty iterator.
        """
        reads = Reads()
        self.assertEqual([], list(reads))

    def testNoReadsLength(self):
        """
        A Reads instance with no reads should have a length of zero.
        """
        reads = Reads()
        self.assertEqual(0, len(reads))

    def testManuallyAddedReads(self):
        """
        A Reads instance with reads added manually must be able to be listed.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        self.assertEqual([read1, read2], list(reads))

    def testManuallyAddedReadsLength(self):
        """
        A Reads instance with reads added manually must have the correct
        length.
        """
        reads = Reads()
        reads.add(Read('id1', 'AT'))
        reads.add(Read('id2', 'AC'))
        self.assertEqual(2, len(reads))

    def testSubclass(self):
        """
        A Reads subclass with an iter method should result in an instance
        with a correct iterator.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')

        class FastaReads(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = FastaReads()
        self.assertEqual([read1, read2], list(reads))

    def testSubclassLength(self):
        """
        A Reads subclass with an iter method should result in an instance
        with a correct length.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')

        class FastaReads(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = FastaReads()
        self.assertEqual(2, len(list(reads)))

    def testRepeatedIter(self):
        """
        A Reads subclass with an iter method should be able to be listed
        more than once.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')

        class FastaReads(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = FastaReads()
        self.assertEqual([read1, read2], list(reads))
        self.assertEqual([read1, read2], list(reads))

    def testLengthAfterRepeatedIter(self):
        """
        A Reads subclass with an iter method should be able to be iterated
        more than once, and following each its length should be correct.
        """
        class FastaReads(Reads):
            def iter(self):
                yield Read('id1', 'AT')
                yield Read('id2', 'AC')

        reads = FastaReads()
        list(reads)
        self.assertEqual(2, len(reads))
        list(reads)
        self.assertEqual(2, len(reads))

    def testSubclassWithAdditionalReads(self):
        """
        A Reads subclass with an iter method that is then added to manually
        should result in an instance with a correct iterator.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        read3 = Read('id3', 'AC')

        class FastaReads(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = FastaReads()
        reads.add(read3)
        self.assertEqual([read1, read2, read3], list(reads))

    def testSaveWithUnknownFormat(self):
        """
        A Reads instance must raise ValueError if asked to save in an unknown
        format.
        """
        reads = Reads()
        read1 = Read('id1', 'AT', '!!')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        error = "Save format must be either 'fasta' or 'fastq'\\."
        self.assertRaisesRegexp(ValueError, error, reads.save, 'file', 'xxx')

    def testSaveFASTAIsDefault(self):
        """
        A Reads instance must save in FASTA format by default.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            reads.save('filename')
        handle = mockOpener()
        self.assertEqual([call('>id1\nAT\n'), call('>id2\nAC\n')],
                         handle.write.call_args_list)

    def testSaveAsFASTA(self):
        """
        A Reads instance must be able to save in FASTA format.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            reads.save('filename', 'fasta')
        handle = mockOpener()
        self.assertEqual([call('>id1\nAT\n'), call('>id2\nAC\n')],
                         handle.write.call_args_list)

    def testSaveWithUppercaseFormat(self):
        """
        A Reads instance must save correctly when the format string is
        given in upper case.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            reads.save('filename', 'FASTA')
        handle = mockOpener()
        self.assertEqual([call('>id1\nAT\n'), call('>id2\nAC\n')],
                         handle.write.call_args_list)

    def testSaveAsFASTQ(self):
        """
        A Reads instance must be able to save in FASTQ format.
        """
        reads = Reads()
        read1 = Read('id1', 'AT', '!!')
        read2 = Read('id2', 'AC', '@@')
        reads.add(read1)
        reads.add(read2)
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            reads.save('filename', 'fastq')
        handle = mockOpener()
        self.assertEqual(
            [call('@id1\nAT\n+id1\n!!\n'), call('@id2\nAC\n+id2\n@@\n')],
            handle.write.call_args_list)

    def testSaveAsFASTQFailsOnReadWithNoQuality(self):
        """
        A Reads instance must raise a ValueError if asked to save in FASTQ
        format and there is a read with no quality present.
        """
        reads = Reads()
        read1 = Read('id1', 'AT', '!!')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        error = "Read 'id2' has no quality information"
        self.assertRaisesRegexp(ValueError, error, reads.save, 'file', 'fastq')

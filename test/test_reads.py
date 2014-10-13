from unittest import TestCase
from mock import patch, call
from cStringIO import StringIO

from mocking import mockOpen
from dark.reads import (Read, TranslatedRead, Reads, PropertiesRead,
                        DNARNARead, AARead)


class TestRead(TestCase):
    """
    Test the Read class.
    """
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

    def testLength(self):
        """
        len() must return the length of a read's sequence.
        """
        read = Read('id', 'ACGT', '!!!!')
        self.assertEqual(4, len(read))

    def testToUnknownFormat(self):
        """
        toString must raise a ValueError if asked to convert to an unknown
        format.
        """
        read = Read('id', 'ACGT', '!!!!')
        error = "Format must be either 'fasta' or 'fastq'\\."
        self.assertRaisesRegexp(ValueError, error, read.toString, 'unknown')

    def testToFASTA(self):
        """
        toString must return correct FASTA.
        """
        read = Read('id', 'ACGT')
        self.assertEqual('>id\nACGT\n', read.toString('fasta'))

    def testToFASTAWithQuality(self):
        """
        toString must return correct FASTA, including when a read has quality
        information (which is not present in FASTA).
        """
        read = Read('id', 'ACGT', '!!!!')
        self.assertEqual('>id\nACGT\n', read.toString('fasta'))

    def testToFASTQWithNoQuality(self):
        """
        toString must raise a ValueError if asked to convert to FASTQ but the
        read has no quality.
        """
        read = Read('id', 'ACGT')
        error = "Read 'id' has no quality information"
        self.assertRaisesRegexp(ValueError, error, read.toString, 'fastq')

    def testToFASTQ(self):
        """
        toString must return correct FASTA.
        """
        read = Read('id', 'ACGT', '!@#$')
        self.assertEqual('@id\nACGT\n+id\n!@#$\n', read.toString('fastq'))

    def testEqualityWithDifferingIds(self):
        """
        If two Read instances have different ids, they must not be
        considered equal.
        """
        self.assertNotEqual(Read('id1', 'AC'), Read('id2', 'AC'))

    def testEqualityWithDifferingSequences(self):
        """
        If two Read instances have different sequences, they must not be
        considered equal.
        """
        self.assertNotEqual(Read('id1', 'AA'), Read('id1', 'CC'))

    def testEqualityWithDifferingQuality(self):
        """
        If two Read instances have different quality, they must not be
        considered equal.
        """
        self.assertNotEqual(Read('id1', 'AC', 'qq'), Read('id1', 'AC', 'rr'))

    def testEqualityWithOneOmittedQuality(self):
        """
        If two Read instances have different quality (one is omitted), they
        must not be considered equal.
        """
        self.assertNotEqual(Read('id1', 'AC'), Read('id1', 'AC', 'rr'))

    def testEqualityWithNoQuality(self):
        """
        If two Read instances have the same id and sequence and neither has a
        quality, they must be considered equal.
        """
        self.assertEqual(Read('id1', 'AC'), Read('id1', 'AC'))

    def testEquality(self):
        """
        If two Read instances have the same id, sequence, and quality, they
        must be considered equal.
        """
        self.assertEqual(Read('id1', 'AC', 'qq'), Read('id1', 'AC', 'qq'))


class TestDNARNARead(TestCase):
    """
    Tests for the DNARNARead class.
    """
    def testReverseComplementReversesQuality(self):
        """
        The reverseComplement function must return a reversed quality string.
        """
        read = DNARNARead('id', 'atcg', quality='!@#$')
        self.assertEqual('$#@!', read.reverseComplement().quality)

    def testReverseComplementDNA(self):
        """
        The reverseComplement function must work for DNA
        """
        read = DNARNARead('id', 'atcg', quality='!@#$')
        self.assertEqual('CGAT', read.reverseComplement().sequence)

    def testReverseComplementAmbiguousDNA(self):
        """
        The reverseComplement function must work for DNA that includes
        ambiguous bases.
        """
        read = DNARNARead('id', 'atcgmrwsvhxn')
        self.assertEqual('NXDBSWYKCGAT', read.reverseComplement().sequence)

    def testReverseComplementRNA(self):
        """
        The reverseComplement function must work for RNA
        """
        read = DNARNARead('id', 'aucg')
        self.assertEqual('CGAU', read.reverseComplement().sequence)

    def testReverseComplementAmbiguousRNA(self):
        """
        The reverseComplement function must work for RNA that includes
        ambiguous bases.
        """
        read = DNARNARead('id', 'aucgmrwsykvhxn')
        self.assertEqual('NXDBMRSWYKCGAU', read.reverseComplement().sequence)

    def testTranslationsOfEmptySequence(self):
        """
        The translations function must correctly return all six (empty)
        translations of the empty sequence.
        """
        read = DNARNARead('id', '')
        self.assertEqual(
            [
                TranslatedRead(read, '', 0, False),
                TranslatedRead(read, '', 1, False),
                TranslatedRead(read, '', 2, False),
                TranslatedRead(read, '', 0, True),
                TranslatedRead(read, '', 1, True),
                TranslatedRead(read, '', 2, True)
            ],
            list(read.translations()))

    def testTranslationsOfOneBaseSequence(self):
        """
        The translations function must correctly return all six translations
        of a sequence with just one base.
        """
        read = DNARNARead('id', 'a')
        self.assertEqual(
            [
                TranslatedRead(read, 'X', 0, False),
                TranslatedRead(read, '', 1, False),
                TranslatedRead(read, '', 2, False),
                TranslatedRead(read, 'X', 0, True),
                TranslatedRead(read, '', 1, True),
                TranslatedRead(read, '', 2, True)
            ],
            list(read.translations()))

    def testTranslationsOfTwoBaseSequence(self):
        """
        The translations function must correctly return all six translations
        of a sequence with just two bases.
        """
        read = DNARNARead('id', 'ag')
        self.assertEqual(
            [
                TranslatedRead(read, 'X', 0, False),
                TranslatedRead(read, 'X', 1, False),
                TranslatedRead(read, '', 2, False),
                TranslatedRead(read, 'L', 0, True),
                TranslatedRead(read, 'X', 1, True),
                TranslatedRead(read, '', 2, True)
            ],
            list(read.translations()))

    def testTranslationOfStopCodonTAA(self):
        """
        The translations function must correctly translate the TAA stop codon.
        """
        read = DNARNARead('id', 'taa')
        self.assertEqual(
            TranslatedRead(read, '*', 0, False),
            read.translations().next())

    def testTranslationOfStopCodonTAG(self):
        """
        The translations function must correctly translate the TAG stop codon.
        """
        read = DNARNARead('id', 'tag')
        self.assertEqual(
            TranslatedRead(read, '*', 0, False),
            read.translations().next())

    def testTranslationOfStopCodonTGA(self):
        """
        The translations function must correctly translate the TGA stop codon.
        """
        read = DNARNARead('id', 'tga')
        self.assertEqual(
            TranslatedRead(read, '*', 0, False),
            read.translations().next())

    def testTranslationOfMultipleStopCodons(self):
        """
        The translations function must correctly translate multiple stop codons
        in a sequence.
        """
        read = DNARNARead('id', 'taatagtga')
        self.assertEqual(
            TranslatedRead(read, '***', 0, False),
            read.translations().next())

    def testTranslationOfStartCodonATG(self):
        """
        The translations function must correctly translate the ATG start codon
        to a methionine (M).
        """
        read = DNARNARead('id', 'atg')
        self.assertEqual(
            TranslatedRead(read, 'M', 0, False),
            read.translations().next())

    def testTranslations(self):
        """
        The translations function must correctly return all six translations.
        """
        read = DNARNARead('id', 'accgtcagg')
        self.assertEqual(
            [
                TranslatedRead(read, 'TVR', 0, False),
                TranslatedRead(read, 'PSG', 1, False),
                TranslatedRead(read, 'RQX', 2, False),
                TranslatedRead(read, 'PDG', 0, True),
                TranslatedRead(read, 'LTV', 1, True),
                TranslatedRead(read, '*RX', 2, True)
            ],
            list(read.translations()))


class TestAARead(TestCase):
    """
    Tests for the AARead class.
    """
    def testAAToPropertiesCorrectTranslation(self):
        """
        An AA sequence must be correctly translated to properties.
        """
        read = AARead('id', 'ADADR*')
        result = read.aaToProperties()
        self.assertEqual(PropertiesRead(read,
                         [193, 3202, 193, 3202, 2562, 4096]), result)


class TestPropertiesRead(TestCase):
    """
    Test the PropertiesRead class.
    """
    def testAllAttributes(self):
        """
        A PropertiesRead instance must have the expected attributes.
        """
        read = AARead('id', 'ADADR')
        properties = PropertiesRead(read, [193, 3202, 193, 3202, 2562])
        self.assertEqual([193, 3202, 193, 3202, 2562], properties.sequence)
        self.assertEqual(read, properties.originalRead)

    def testRightSequence(self):
        """
        A PropertiesRead must have the expected sequence.
        """
        read = AARead('id', 'ADADR')
        properties = PropertiesRead(read, [193, 3202, 193, 3202, 2562])
        self.assertEqual([193, 3202, 193, 3202, 2562], properties.sequence)

    def testRightId(self):
        """
        A PropertiesRead must have the right id.
        """
        read = AARead('id', 'ADADR')
        properties = PropertiesRead(read, 'HAHAB')
        self.assertEqual('id-properties', properties.id)


class TestTranslatedRead(TestCase):
    """
    Test the TranslatedRead class.
    """

    def testExpectedAttributes(self):
        """
        A TranslatedRead instance must have the expected attributes.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 0)
        self.assertEqual('IRDS', translated.sequence)
        self.assertEqual(0, translated.frame)
        self.assertIs(read, translated.originalRead)

    def testSequence(self):
        """
        A TranslatedRead instance must have the expected sequence.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 0)
        self.assertEqual('IRDS', translated.sequence)

    def testOutOfRangeFrame(self):
        """
        A TranslatedRead instance must raise a ValueError if the passed frame
        is not 0, 1, or 2.
        """
        read = Read('id', 'atcgatcgatcg')
        error = 'Frame must be 0, 1, or 2'
        self.assertRaisesRegexp(ValueError, error, TranslatedRead, read,
                                'IRDS', 3)

    def testExpectedFrame(self):
        """
        A TranslatedRead instance must have the expected frame.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 2)
        self.assertEqual(2, translated.frame)

    def testReverseComplemented(self):
        """
        A TranslatedRead instance must have the expected reversedComplemented
        value.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 0)
        self.assertFalse(translated.reverseComplemented)
        translated = TranslatedRead(read, 'IRDS', 0, reverseComplemented=True)
        self.assertTrue(translated.reverseComplemented)

    def testExpectedOriginalRead(self):
        """
        A TranslatedRead instance must store the original read.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 0)
        self.assertIs(read, translated.originalRead)

    def testId(self):
        """
        A TranslatedRead instance must put the the frame information into its
        read id.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 0)
        self.assertEqual('id-frame0', translated.id)

    def testIdReverseComplemented(self):
        """
        A TranslatedRead instance must put the the frame information into its
        read id when the original read was reverse complemented.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 1, True)
        self.assertEqual('id-frame1rc', translated.id)

    def testMaximumORFLengthNoStops(self):
        """
        The maximumORFLength function must return the correct value when
        there are no stop codons in a translated read.
        """
        read = Read('id', 'atcgatcgatcg')
        translated = TranslatedRead(read, 'IRDS', 0)
        self.assertEqual(4, translated.maximumORFLength())

    def testMaximumORFLength(self):
        """
        The maximumORFLength function must return the correct value.
        """
        read = Read('id', 'acctaggttgtttag')
        translated = TranslatedRead(read, 'T*VV*', 0)
        self.assertEqual(2, translated.maximumORFLength())


class TestReads(TestCase):
    """
    Test the Reads class.
    """

    def testNoReads(self):
        """
        A Reads instance with no reads must return an empty iterator.
        """
        reads = Reads()
        self.assertEqual([], list(reads))

    def testNoReadsLength(self):
        """
        A Reads instance with no reads must have a length of zero.
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
        A Reads subclass with an iter method must result in an instance
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
        A Reads subclass with an iter method must result in an instance
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
        A Reads subclass with an iter method must be able to be listed
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
        A Reads subclass with an iter method must be able to be iterated
        more than once, and following each its length must be correct.
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
        must result in an instance with a correct iterator.
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
        error = "Format must be either 'fasta' or 'fastq'\\."
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

    def testSaveToFileDescriptor(self):
        """
        A Reads instance must save to a file-like object if not passed a string
        filename.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        fp = StringIO()
        reads.save(fp)
        self.assertEqual('>id1\nAT\n>id2\nAC\n', fp.getvalue())

    def testFilterOnMinLength(self):
        """
        Filtering on minimal length must work.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ACG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(minLength=4)
        self.assertEqual([read1], list(result))

    def testFilterOnMaxLength(self):
        """
        Filtering on maximal length must work.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ACG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(maxLength=3)
        self.assertEqual([read2], list(result))

    def testFilterOnLengthNothingMatches(self):
        """
        When filtering on length, no reads must be returned if none of them
        satisfy the length requirements.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ACG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(minLength=10, maxLength=15)
        self.assertEqual([], list(result))

    def testFilterOnLengthEverythingMatches(self):
        """
        When filtering on length, all reads must be returned if they all
        satisfy the length requirements.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ACG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(minLength=2, maxLength=5)
        self.assertEqual([read1, read2], list(result))

    def testFilterWithMinLengthEqualToMaxLength(self):
        """
        When filtering on length, a read must be returned if its length
        equals a passed minimum and maximum length.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ACG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(minLength=4, maxLength=4)
        self.assertEqual([read1], list(result))

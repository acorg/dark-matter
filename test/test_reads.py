from unittest import TestCase
from unittest.mock import patch, call
from io import StringIO
from os import stat
import builtins

from .mocking import mockOpen
from dark.reads import (
    Read, TranslatedRead, Reads, DNARead, RNARead, AARead, AAReadORF,
    AAReadWithX)
from dark.aa import (
    BASIC_POSITIVE, HYDROPHOBIC, HYDROPHILIC, NEGATIVE, NONE, POLAR, SMALL,
    TINY)
from dark.hsp import HSP


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
        with self.assertRaisesRegex(ValueError, error):
            Read('id', 'ACGT', '!!!')

    def testNoQuality(self):
        """
        If no quality information is given, the read's 'quality' attribute must
        be None.
        """
        read = Read('id', 'ACGT')
        self.assertIs(None, read.quality)

    def testCasePreservation(self):
        """
        The sequence passed to Read must not have its case converted.
        """
        read = Read('id', 'aCGt')
        self.assertEqual('aCGt', read.sequence)

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
        self.assertRaisesRegex(ValueError, error, read.toString, 'unknown')

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
        self.assertRaisesRegex(ValueError, error, read.toString, 'fastq')

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

    def testHashDiffersIfIdDiffers(self):
        """
        The __hash__ value for two reads must differ if their ids differ.
        """
        self.assertNotEqual(hash(Read('id1', 'AA')),
                            hash(Read('id2', 'AA')))

    def testHashDiffersIfSequenceDiffers(self):
        """
        The __hash__ value for two reads must differ if their sequences
        differ.
        """
        self.assertNotEqual(hash(Read('id', 'AAT')),
                            hash(Read('id', 'AAG')))

    def testHashDiffersIfQualityDiffers(self):
        """
        The __hash__ value for two reads must differ if their quality strings
        differ.
        """
        self.assertNotEqual(hash(Read('id', 'AA', '!!')),
                            hash(Read('id', 'AA', '++')))

    def testHashIdenticalNoQuality(self):
        """
        The __hash__ value for two identical reads (with no quality strings)
        must be identical.
        """
        self.assertEqual(hash(Read('id', 'AA')),
                         hash(Read('id', 'AA')))

    def testHashIdenticalWithQuality(self):
        """
        The __hash__ value for two identical reads (with quality strings) must
        be identical.
        """
        self.assertEqual(hash(Read('id', 'AA', '!!')),
                         hash(Read('id', 'AA', '!!')))

    def testLowComplexityFractionEmptySequence(self):
        """
        A read with an empty sequence must return a zero result from its
        lowComplexityFraction method.
        """
        read = Read('id', '')
        self.assertEqual(0.0, read.lowComplexityFraction())

    def testLowComplexityFractionZero(self):
        """
        A read with no low-complexity bases must return a zero result from its
        lowComplexityFraction method.
        """
        read = Read('id', 'ACGT')
        self.assertEqual(0.0, read.lowComplexityFraction())

    def testLowComplexityFractionOne(self):
        """
        A read with all low-complexity bases must return a one result from its
        lowComplexityFraction method.
        """
        read = Read('id', 'acgt')
        self.assertEqual(1.0, read.lowComplexityFraction())

    def testLowComplexityFraction(self):
        """
        A read with all low-complexity bases must return the correct result
        from its lowComplexityFraction method.
        """
        read = Read('id', 'aCGT')
        self.assertEqual(0.25, read.lowComplexityFraction())

    def testWalkHSPExactMatch(self):
        """
        If the HSP specifies that the entire read matches the subject exactly,
        walkHSP must return the correct results.

        Subject:     ACGT
        Read:        ACGT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=0, readEnd=4, readStartInSubject=0,
                  readEndInSubject=4, subjectStart=0, subjectEnd=4,
                  readMatchedSequence='ACGT', subjectMatchedSequence='ACGT')
        self.assertEqual([(0, 'A', True),
                          (1, 'C', True),
                          (2, 'G', True),
                          (3, 'T', True)],
                         list(read.walkHSP(hsp)))

    def testWalkHSPExactMatchWithGap(self):
        """
        If the HSP specifies that the entire read matches the subject exactly,
        including a gap, walkHSP must return the correct results.

        Subject:     ACGT
        Read:        A-GT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=0, readEnd=4, readStartInSubject=0,
                  readEndInSubject=4, subjectStart=0, subjectEnd=4,
                  readMatchedSequence='A-GT', subjectMatchedSequence='ACGT')
        self.assertEqual([(0, 'A', True),
                          (1, '-', True),
                          (2, 'G', True),
                          (3, 'T', True)],
                         list(read.walkHSP(hsp)))

    def testWalkHSPLeftOverhangingMatch(self):
        """
        If the HSP specifies that the entire read matches the subject, and
        also extends to the left of the subject, walkHSP must return the
        correct results.

        Subject:       GT.....
        Read:        ACGT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=2, readEnd=4, readStartInSubject=-2,
                  readEndInSubject=2, subjectStart=0, subjectEnd=2,
                  readMatchedSequence='GT', subjectMatchedSequence='GT')
        self.assertEqual([(-2, 'A', False),
                          (-1, 'C', False),
                          (0, 'G', True),
                          (1, 'T', True)],
                         list(read.walkHSP(hsp)))

    def testWalkHSPRightOverhangingMatch(self):
        """
        If the HSP specifies that the entire read matches the subject, and
        also extends to the right of the subject, walkHSP must return the
        correct results.

        Subject:       AC
        Read:          ACGT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=0, readEnd=2, readStartInSubject=10,
                  readEndInSubject=14, subjectStart=10, subjectEnd=12,
                  readMatchedSequence='AC', subjectMatchedSequence='AC')
        self.assertEqual([(10, 'A', True),
                          (11, 'C', True),
                          (12, 'G', False),
                          (13, 'T', False)],
                         list(read.walkHSP(hsp)))

    def testWalkHSPLeftAndRightOverhangingMatch(self):
        """
        If the HSP specifies that the read matches the entire subject, and
        also extends to both the left and right of the subject, walkHSP must
        return the correct results.

        Subject:        CG
        Read:          ACGT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=1, readEnd=3, readStartInSubject=10,
                  readEndInSubject=14, subjectStart=11, subjectEnd=13,
                  readMatchedSequence='CG', subjectMatchedSequence='CG')
        self.assertEqual([(10, 'A', False),
                          (11, 'C', True),
                          (12, 'G', True),
                          (13, 'T', False)],
                         list(read.walkHSP(hsp)))

    def testcheckAlphabetwithReadMustBePermissive(self):
        """
        The checkAlphabet function must be permissive if a dark.Read is
        passed.
        """
        read = Read('id', 'ARSTGATGCASASASASASAS')
        self.assertEqual(set('ACGSRT'), read.checkAlphabet())

    def testcheckAlphabetAAReadMatchingReturnTrue(self):
        """
        If an AA read with an AARead readClass is passed in, the checkAlphabet
        function must return the alphabet of the sequence.
        """
        read = AARead('id', 'ARSTGATGCASASASASASAS')
        self.assertEqual(set('ACGSRT'), read.checkAlphabet())

    def testcheckAlphabetDNAReadMatchingReturnTrue(self):
        """
        If a DNA read with a DNARead readClass is passed in, the checkAlphabet
        function must return the alphabet of the sequence.
        """
        read = DNARead('id', 'AAATTAACGGGCCTAGG')
        self.assertEqual(set('ACTG'), read.checkAlphabet())

    def testcheckAlphabetAAReadNotMatchingRaise(self):
        """
        If an AA read with a DNARead readClass is passed in, the checkAlphabet
        function must raise an IndexError.
        """
        read = AARead('id', 'AAATTAACGGGCCTAGG')
        error = "It looks like a DNA sequence has been passed to AARead()."
        self.assertRaisesRegex(ValueError, error, read.checkAlphabet)

    def testcheckAlphabetDNAReadNotMatchingRaise(self):
        """
        If a DNA read with an AARead readClass is passed in, the checkAlphabet
        function must raise an IndexError.
        """
        read = DNARead('id', 'ARSTGATGCASASASASASAS')
        error = ("Read alphabet \('ACGRST'\) is not a subset of expected "
                 "alphabet \('ACGT'\) for read class DNARead.")
        self.assertRaisesRegex(ValueError, error, read.checkAlphabet)


class TestDNARead(TestCase):
    """
    Tests for the DNARead class.
    """
    def testReverseComplementReversesQuality(self):
        """
        The reverseComplement function must return a reversed quality string.
        """
        read = DNARead('id', 'atcg', quality='!@#$')
        self.assertEqual('$#@!', read.reverseComplement().quality)

    def testReverseComplement(self):
        """
        The reverseComplement function must work.
        """
        read = DNARead('id', 'ATCG', quality='!@#$')
        self.assertEqual('CGAT', read.reverseComplement().sequence)

    def testReverseComplementAmbiguous(self):
        """
        The reverseComplement function must work for a sequence that includes
        ambiguous bases.
        """
        read = DNARead('id', 'ATCGMRWSVHXN')
        self.assertEqual('NXDBSWYKCGAT', read.reverseComplement().sequence)

    def testTranslationsOfEmptySequence(self):
        """
        The translations function must correctly return all six (empty)
        translations of the empty sequence.
        """
        read = DNARead('id', '')
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
        read = DNARead('id', 'a')
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
        read = DNARead('id', 'AG')
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

    def testTranslationOfStopCodonTAG(self):
        """
        The translations function must correctly translate the TAG stop codon.
        """
        read = DNARead('id', 'tag')
        self.assertEqual(
            TranslatedRead(read, '*', 0, False),
            next(read.translations()))

    def testTranslationOfStopCodonTGA(self):
        """
        The translations function must correctly translate the TGA stop codon.
        """
        read = DNARead('id', 'tga')
        self.assertEqual(
            TranslatedRead(read, '*', 0, False),
            next(read.translations()))

    def testTranslationOfMultipleStopCodons(self):
        """
        The translations function must correctly translate multiple stop codons
        in a sequence.
        """
        read = DNARead('id', 'taatagtga')
        self.assertEqual(
            TranslatedRead(read, '***', 0, False),
            next(read.translations()))

    def testTranslationOfStartCodonATG(self):
        """
        The translations function must correctly translate the ATG start codon
        to a methionine (M).
        """
        read = DNARead('id', 'atg')
        self.assertEqual(
            TranslatedRead(read, 'M', 0, False),
            next(read.translations()))

    def testTranslations(self):
        """
        The translations function must correctly return all six translations.
        """
        read = DNARead('id', 'ACCGTCAGG')
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


class TestRNARead(TestCase):
    """
    Tests for the RNARead class.
    """
    def testReverseComplement(self):
        """
        The reverseComplement function must work.
        """
        read = RNARead('id', 'AUCG')
        self.assertEqual('CGAU', read.reverseComplement().sequence)

    def testReverseComplementAmbiguous(self):
        """
        The reverseComplement function must work for a sequence that includes
        ambiguous bases.
        """
        read = RNARead('id', 'AUCGMRWSYKVHXN')
        self.assertEqual('NXDBMRSWYKCGAU', read.reverseComplement().sequence)

    def testTranslationOfStopCodonUAA(self):
        """
        The translations function must correctly translate the UAA stop codon.
        """
        read = RNARead('id', 'UAA')
        self.assertEqual(
            TranslatedRead(read, '*', 0, False),
            next(read.translations()))


class TestAARead(TestCase):
    """
    Tests for the AARead class.
    """
    def testPropertiesCorrectTranslation(self):
        """
        The properties function must work correctly.
        """
        read = AARead('id', 'ADR*')
        properties = read.properties()
        self.assertEqual(
            [
                HYDROPHOBIC | SMALL | TINY,
                HYDROPHILIC | SMALL | POLAR | NEGATIVE,
                HYDROPHILIC | POLAR | BASIC_POSITIVE,
                NONE
            ],
            list(properties)
        )

    def testPropertyDetailsCorrectTranslation(self):
        """
        The propertyDetails function must return the right property details
        sequence.
        """
        read = AARead('id', 'ADR*')
        properties = read.propertyDetails()
        self.assertEqual(
            [
                {
                    'polarity': -0.20987654321,
                    'aliphaticity': 0.305785123967,
                    'volume': -0.664670658683,
                    'polar_req': -0.463414634146,
                    'hydropathy': 0.4,
                    'iep': -0.191489361702,
                    'hydroxyethilation': -0.265160523187,
                    'aromaticity': -0.550128534704,
                    'hydrogenation': 0.8973042362,
                    'composition': -1.0
                },
                {
                    'polarity': 1.0,
                    'aliphaticity': -0.818181818182,
                    'volume': -0.389221556886,
                    'polar_req': 1.0,
                    'hydropathy': -0.777777777778,
                    'iep': -1.0,
                    'hydroxyethilation': -0.348394768133,
                    'aromaticity': -1.0,
                    'hydrogenation': -0.90243902439,
                    'composition': 0.00363636363636
                },
                {
                    'polarity': 0.382716049383,
                    'aliphaticity': -0.157024793388,
                    'volume': 0.449101796407,
                    'polar_req': 0.0487804878049,
                    'hydropathy': -1.0,
                    'iep': 1.0,
                    'hydroxyethilation': -0.51486325802,
                    'aromaticity': -0.0642673521851,
                    'hydrogenation': -0.401797175866,
                    'composition': -0.527272727273
                },
                NONE],
            list(properties)
        )

    def testORFsEmptySequence(self):
        """
        An AA read of length zero must not have any ORFs.
        """
        read = AARead('id', '')
        orfs = list(read.ORFs())
        self.assertEqual(0, len(orfs))

    def testORFsEmptySequenceWithStartStop(self):
        """
        An AA read with just a start and stop codon must not have any ORFs.
        """
        read = AARead('id', 'M*')
        orfs = list(read.ORFs())
        self.assertEqual(0, len(orfs))

    def testORFsEmptySequenceWithStart(self):
        """
        An AA read with just a start codon must not have any ORFs.
        """
        read = AARead('id', 'M')
        orfs = list(read.ORFs())
        self.assertEqual(0, len(orfs))

    def testORFsWithOneStopCodon(self):
        """
        An AA read of a single stop codon must not have any ORFs.
        """
        read = AARead('id', '*')
        orfs = list(read.ORFs())
        self.assertEqual(0, len(orfs))

    def testORFsWithTwoStopCodons(self):
        """
        An AA read of two stop codons must not have any ORFs.
        """
        read = AARead('id', '**')
        orfs = list(read.ORFs())
        self.assertEqual(0, len(orfs))

    def testORFsWithJustStartsAndStops(self):
        """
        An AA read of only start and stop codons must not have any ORFs.
        """
        read = AARead('id', '**MM*M**MMM*')
        orfs = list(read.ORFs())
        self.assertEqual(0, len(orfs))

    def testOpenOpenORF(self):
        """
        An AA read that contains no start or stop codons should result in
        just one AAReadORF when its ORFs method is called. The ORF must have
        the correct start/stop offsets and its left and right side must be
        marked as open.
        """
        read = AARead('id', 'ADRADR')
        orfs = list(read.ORFs())
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(6, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:6)')

    def testOpenCloseORF(self):
        """
        An AA read that contains no start codon but one trailing stop
        codon should result in just one AAReadORF when its ORFs method is
        called. The ORF must have the correct start/stop offsets and its
        left and right sides must be marked as open and closed,
        respectively.
        """
        read = AARead('id', 'ADRADR*')
        orfs = list(read.ORFs())
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(6, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:6]')

    def testOpenCloseORFWithMultipleStops(self):
        """
        An AA read that contains no start codon but multiple trailing
        stop codons should result in just one AAReadORF when its ORFs
        method is called. The ORF must have the correct start/stop offsets
        and its left and right sides must be marked as open and closed,
        respectively.
        """
        read = AARead('id', 'ADRADR***')
        orfs = list(read.ORFs())
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(6, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:6]')

    def testCloseOpenORF(self):
        """
        An AA read that contains a start but no stop codon should result in
        just one AAReadORF when its ORFs method is called. The ORF must have
        the correct start/stop offsets and its left and right sides must be
        marked as closed and open, respectively.
        """
        read = AARead('id', 'MADRADR')
        orfs = list(read.ORFs())
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(1, orf.start)
        self.assertEqual(7, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[1:7)')

    def testCloseOpenORFWithMultipleStarts(self):
        """
        An AA read that contains multiple initial start codons but no
        stop codon should result in just one AAReadORF when its ORFs method
        is called. The ORF must have the correct start/stop offsets and its
        left and right sides must be marked as closed and open,
        respectively.
        """
        read = AARead('id', 'MMMADRADR')
        orfs = list(read.ORFs())
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(3, orf.start)
        self.assertEqual(9, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[3:9)')

    def testCloseCloseORF(self):
        """
        An AA read that contains a start and a stop codon should result in
        just one AAReadORF when its ORFs method is called. The ORF must have
        the correct start/stop offsets and its left and right sides must be
        both marked as closed.
        """
        read = AARead('id', 'MADRADR*')
        orfs = list(read.ORFs())
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(1, orf.start)
        self.assertEqual(7, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[1:7]')

    def testOpenCloseThenCloseOpenORF(self):
        """
        An AA read that contains an ORF that is left open and right closed
        followed by an ORF that is left closed and right open must have the
        ORFs detected correctly when its ORFs method is called.
        """
        read = AARead('id', 'ADR*MRRR')
        orfs = list(read.ORFs())
        self.assertEqual(2, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(3, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:3]')

        orf = orfs[1]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(5, orf.start)
        self.assertEqual(8, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[5:8)')

    def testCloseCloseThenCloseOpenORF(self):
        """
        An AA read that contains an ORF that is left and right closed
        followed by an ORF that is left closed and right open must have the
        ORFs detected correctly when its ORFs method is called.
        """
        read = AARead('id', 'MADR*MRRR')
        orfs = list(read.ORFs())
        self.assertEqual(2, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(1, orf.start)
        self.assertEqual(4, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[1:4]')

        orf = orfs[1]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(6, orf.start)
        self.assertEqual(9, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[6:9)')

    def testCloseCloseThenCloseCloseORF(self):
        """
        An AA read that contains two ORFs that are both left and right closed
        must have the ORFs detected correctly when its ORFs method is called.
        """
        read = AARead('id', 'MADR*MRRR*')
        orfs = list(read.ORFs())
        self.assertEqual(2, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(1, orf.start)
        self.assertEqual(4, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[1:4]')

        orf = orfs[1]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(6, orf.start)
        self.assertEqual(9, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[6:9]')

    def testOpenCloseThenCloseCloseThenCloseOpenORF(self):
        """
        An AA read that contains an ORF that is left open and right closed
        followed by an internal ORF, followed by an ORF that is left closed
        and right open must have the ORFs detected correctly when its ORFs
        method is called.
        """
        read = AARead('id', 'ADR*MAAA*MRRR')
        orfs = list(read.ORFs())
        self.assertEqual(3, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(3, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:3]')

        orf = orfs[1]
        self.assertEqual('AAA', orf.sequence)
        self.assertEqual(5, orf.start)
        self.assertEqual(8, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[5:8]')

        orf = orfs[2]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(10, orf.start)
        self.assertEqual(13, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[10:13)')

    def testCloseCloseThenCloseCloseThenCloseOpenORF(self):
        """
        An AA read that contains an ORF that is left and right closed
        followed by an internal ORF, followed by an ORF that is left closed
        and right open must have the ORFs detected correctly when its ORFs
        method is called.
        """
        read = AARead('id', 'MADR*MAAA*MRRR')
        orfs = list(read.ORFs())
        self.assertEqual(3, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(1, orf.start)
        self.assertEqual(4, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[1:4]')

        orf = orfs[1]
        self.assertEqual('AAA', orf.sequence)
        self.assertEqual(6, orf.start)
        self.assertEqual(9, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[6:9]')

        orf = orfs[2]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(11, orf.start)
        self.assertEqual(14, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[11:14)')

    def testCloseCloseThenCloseCloseThenCloseCloseORF(self):
        """
        An AA read that contains an ORF that is left and right closed
        followed by an internal ORF, followed by an ORF that is left and
        right closed must have the ORFs detected correctly when its ORFs
        method is called.
        """
        read = AARead('id', 'MADR*MAAA*MRRR*')
        orfs = list(read.ORFs())
        self.assertEqual(3, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(1, orf.start)
        self.assertEqual(4, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[1:4]')

        orf = orfs[1]
        self.assertEqual('AAA', orf.sequence)
        self.assertEqual(6, orf.start)
        self.assertEqual(9, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[6:9]')

        orf = orfs[2]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(11, orf.start)
        self.assertEqual(14, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[11:14]')

    def testOpenCloseThenCloseCloseThenCloseOpenORFWithJunk(self):
        """
        An AA read that contains an ORF that is left open and right closed
        followed by an internal ORF, followed by an ORF that is left closed
        and right open must have the ORFs detected correctly when its ORFs
        method is called, including when there are intermediate start and
        stop codons.
        """
        read = AARead('id', 'ADR***M*MAAA***MMM*MRRR')
        orfs = list(read.ORFs())
        self.assertEqual(3, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(3, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:3]')

        orf = orfs[1]
        self.assertEqual('AAA', orf.sequence)
        self.assertEqual(9, orf.start)
        self.assertEqual(12, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[9:12]')

        orf = orfs[2]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(20, orf.start)
        self.assertEqual(23, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[20:23)')

    def testCloseCloseThenCloseCloseThenCloseOpenORFWithJunk(self):
        """
        An AA read that contains an ORF that is left and right closed
        followed by an internal ORF, followed by an ORF that is left closed
        and right open must have the ORFs detected correctly when its ORFs
        method is called.
        """
        read = AARead('id', '**MADR***MM**MAAA***M*MRRR')
        orfs = list(read.ORFs())
        self.assertEqual(3, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(3, orf.start)
        self.assertEqual(6, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[3:6]')

        orf = orfs[1]
        self.assertEqual('AAA', orf.sequence)
        self.assertEqual(14, orf.start)
        self.assertEqual(17, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[14:17]')

        orf = orfs[2]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(23, orf.start)
        self.assertEqual(26, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[23:26)')

    def testCloseCloseThenCloseCloseThenCloseCloseORFWithJunk(self):
        """
        An AA read that contains an ORF that is left and right closed
        followed by an internal ORF, followed by an ORF that is left and
        right closed must have the ORFs detected correctly when its ORFs
        method is called.
        """
        read = AARead('id', 'M***MADR***MAAA***MRRR***MM')
        orfs = list(read.ORFs())
        self.assertEqual(3, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(5, orf.start)
        self.assertEqual(8, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[5:8]')

        orf = orfs[1]
        self.assertEqual('AAA', orf.sequence)
        self.assertEqual(12, orf.start)
        self.assertEqual(15, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[12:15]')

        orf = orfs[2]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(19, orf.start)
        self.assertEqual(22, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[19:22]')

    def testNoStartCodon_GithubIssue239(self):
        """
        If there is no start codon in a sequence, it should not be returned
        as an ORF.
        """
        # Example from https://github.com/acorg/dark-matter/issues/239
        read = AARead('id', 'KK*LLILFSCQRWSRKSICVHLTQR*G*')
        orfs = list(read.ORFs())
        self.assertEqual(1, len(orfs))

        orf = orfs[0]
        self.assertEqual('KK', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(2, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:2]')

    def testGOR4(self):
        """
        The gor4 method must return the expected predictions and prediction
        probabilities for an AARead.
        """
        # Note that this test just re-does what the testApoamicyanin test
        # checks in test/test_gor4.py but on an AARead instead of on a
        # raw AA sequence string.
        seq = 'DKATIPSESPFAAAEVADGAIVVDIAKMKYETPELHVKVGDTVTWINREA'
        read = AARead('id', seq)
        result = read.gor4()
        expected = 'CCCCCCCCCCHHHHHHHCCHHHHHHHHHHHCCCCEEEEECCEEEEEEEEC'
        self.assertEqual(expected, result['predictions'])

        # Below are the probabilities we get from running the command-line
        # version of the GOR IV code (with the printout function adjusted
        # to print 7 decimal places).
        expected = (
            (0.0000000, 0.0039858, 0.9960142),
            (0.0000000, 0.0278830, 0.9721170),
            (0.0000001, 0.0272379, 0.9727620),
            (0.0000012, 0.1169491, 0.8830497),
            (0.0001399, 0.0341633, 0.9656968),
            (0.0017102, 0.0208701, 0.9774197),
            (0.0105059, 0.0851670, 0.9043271),
            (0.0555277, 0.0633068, 0.8811655),
            (0.1095295, 0.0324757, 0.8579948),
            (0.4153768, 0.0191681, 0.5654551),
            (0.5536253, 0.1005364, 0.3458383),
            (0.6340716, 0.1182708, 0.2476575),
            (0.7778620, 0.0700642, 0.1520738),
            (0.8236157, 0.0761598, 0.1002245),
            (0.7757733, 0.0764664, 0.1477603),
            (0.7078716, 0.0882713, 0.2038570),
            (0.6208869, 0.0486473, 0.3304658),
            (0.3671274, 0.0182623, 0.6146103),
            (0.4309160, 0.0364063, 0.5326777),
            (0.5071787, 0.1901496, 0.3026716),
            (0.4161928, 0.4818794, 0.1019277),
            (0.4236001, 0.5074177, 0.0689823),
            (0.5267394, 0.4115629, 0.0616978),
            (0.6260716, 0.2516718, 0.1222566),
            (0.6345364, 0.2743220, 0.0911416),
            (0.6292611, 0.2126158, 0.1581231),
            (0.5834816, 0.0938588, 0.3226597),
            (0.5517939, 0.0403667, 0.4078394),
            (0.4721551, 0.0882974, 0.4395475),
            (0.4599132, 0.0910693, 0.4490175),
            (0.1938535, 0.0556542, 0.7504924),
            (0.0518396, 0.0844042, 0.8637562),
            (0.1642072, 0.0529393, 0.7828534),
            (0.2086895, 0.3441093, 0.4472012),
            (0.1385842, 0.6654806, 0.1959352),
            (0.0953126, 0.6765780, 0.2281094),
            (0.0322423, 0.9027448, 0.0650129),
            (0.0787300, 0.6267382, 0.2945318),
            (0.0786274, 0.5764974, 0.3448752),
            (0.0711691, 0.1778753, 0.7509556),
            (0.0765578, 0.1883003, 0.7351419),
            (0.1004260, 0.7001448, 0.1994292),
            (0.0736500, 0.8145193, 0.1118307),
            (0.0321693, 0.8462602, 0.1215706),
            (0.0968932, 0.6217796, 0.2813272),
            (0.0629602, 0.7199594, 0.2170804),
            (0.0186437, 0.5259689, 0.4553874),
            (0.0009734, 0.5355605, 0.4634661),
            (0.0000139, 0.8630118, 0.1369742),
            (0.0000002, 0.9903608, 0.0096390))

        for e, r in zip(expected, result['probabilities']):
            self.assertAlmostEqual(e[0], r[0])
            self.assertAlmostEqual(e[1], r[1])
            self.assertAlmostEqual(e[2], r[2])


class TestAAReadWithX(TestCase):
    """
    Tests for the TestAAReadWithX class.
    """
    def testAlphabet(self):
        """
        The correct alphabet must be used.
        """
        read = AAReadWithX('id', 'ATFDX')
        expected = set('ACDEFGHIKLMNPQRSTVWXY')
        self.assertEqual(expected, read.ALPHABET)

    def testAlphabetChecking(self):
        """
        The alphabet check must work.
        """
        read = AAReadWithX('id', 'ARDGGCFFXEE')
        self.assertEqual(set('ARDGCFFXE'), read.checkAlphabet())


class TestAAReadORF(TestCase):
    """
    Test the AAReadORF class.
    """
    def testSequence(self):
        """
        An AAReadORF instance must have the correct sequence.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 0, 4, True, True)
        self.assertEqual('ADRA', read.sequence)

    def testStart(self):
        """
        An AAReadORF instance must store the correct start offset.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 3, 4, True, True)
        self.assertEqual(3, read.start)
        read = AAReadORF(originalRead, 4, 4, True, True)
        self.assertEqual(4, read.start)

    def testStop(self):
        """
        An AAReadORF instance must store the correct stop offset.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 0, 4, True, True)
        self.assertEqual(4, read.stop)
        read = AAReadORF(originalRead, 0, 6, True, True)
        self.assertEqual(6, read.stop)

    def testOpenLeft(self):
        """
        An AAReadORF instance must store the correct openLeft value.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 0, 4, True, True)
        self.assertTrue(read.openLeft)
        read = AAReadORF(originalRead, 0, 4, False, True)
        self.assertFalse(read.openLeft)

    def testOpenRight(self):
        """
        An AAReadORF instance must store the correct openRight value.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 0, 4, True, True)
        self.assertTrue(read.openRight)
        read = AAReadORF(originalRead, 0, 4, True, False)
        self.assertFalse(read.openRight)

    def testStartGreaterThanStop(self):
        """
        An AAReadORF start offset must not be greater than its stop offset.
        """
        originalRead = AARead('id', 'ADRADR')
        error = 'start offset \(4\) greater than stop offset \(0\)'
        self.assertRaisesRegex(
            ValueError, error, AAReadORF, originalRead, 4, 0, True, True)

    def testStartNegative(self):
        """
        An AAReadORF start offset must not be less than zero.
        """
        originalRead = AARead('id', 'ADRADR')
        error = 'start offset \(-1\) less than zero'
        self.assertRaisesRegex(
            ValueError, error, AAReadORF, originalRead, -1, 6, True, True)

    def testStopGreaterThanOriginalSequenceLength(self):
        """
        An AAReadORF stop offset must not be greater than the length of the
        original sequence.
        """
        originalRead = AARead('id', 'ADRADR')
        error = 'stop offset \(10\) > original read length \(6\)'
        self.assertRaisesRegex(
            ValueError, error, AAReadORF, originalRead, 0, 10, True, True)

    def testOpenOpenId(self):
        """
        An AAReadORF instance must have a correctly annotated (containing the
        sequence offsets) id when the left and right sides of the ORF are open.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 3, 4, True, True)
        self.assertEqual('id-(3:4)', read.id)

    def testOpenClosedId(self):
        """
        An AAReadORF instance must have a correctly annotated (containing the
        sequence offsets) id when the left side of the ORF is open and the
        right is closed.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 3, 4, True, False)
        self.assertEqual('id-(3:4]', read.id)

    def testClosedOpenId(self):
        """
        An AAReadORF instance must have a correctly annotated (containing the
        sequence offsets) id when the left side of the ORF is closed and the
        right is open.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 3, 4, False, True)
        self.assertEqual('id-[3:4)', read.id)

    def testClosedClosedId(self):
        """
        An AAReadORF instance must have a correctly annotated (containing the
        sequence offsets) id when both sides of the ORF are closed.
        """
        originalRead = AARead('id', 'ADRADR')
        read = AAReadORF(originalRead, 3, 4, False, False)
        self.assertEqual('id-[3:4]', read.id)


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
        self.assertRaisesRegex(ValueError, error, TranslatedRead, read,
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
        translated = TranslatedRead(read, 'T*MVV*', 0)
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
        self.assertRaisesRegex(ValueError, error, reads.save, 'file', 'xxx')
        # The output file must not exist following the save() failure.
        error = "No such file or directory: 'file'"
        self.assertRaisesRegex(OSError, error, stat, 'file')

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
        with patch.object(builtins, 'open', mockOpener):
            reads.save('filename')
        handle = mockOpener()
        self.assertEqual([call('>id1\nAT\n'), call('>id2\nAC\n')],
                         handle.write.mock_calls)

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
        with patch.object(builtins, 'open', mockOpener):
            reads.save('filename', 'fasta')
        handle = mockOpener()
        self.assertEqual([call('>id1\nAT\n'), call('>id2\nAC\n')],
                         handle.write.mock_calls)

    def testSaveReturnsReadsInstance(self):
        """
        The save method on a Reads instance must return that instance.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            result = reads.save('filename')
            self.assertIs(reads, result)

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
        with patch.object(builtins, 'open', mockOpener):
            reads.save('filename', 'FASTA')
        handle = mockOpener()
        self.assertEqual([call('>id1\nAT\n'), call('>id2\nAC\n')],
                         handle.write.mock_calls)

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
        with patch.object(builtins, 'open', mockOpener):
            reads.save('filename', 'fastq')
        handle = mockOpener()
        self.assertEqual(
            [call('@id1\nAT\n+id1\n!!\n'), call('@id2\nAC\n+id2\n@@\n')],
            handle.write.mock_calls)

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
        self.assertRaisesRegex(ValueError, error, reads.save, 'file', 'fastq')
        # The output file must not exist following the save() failure.
        error = "No such file or directory: 'file'"
        self.assertRaisesRegex(OSError, error, stat, 'file')

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

    def testFilterNoArgs(self):
        """
        Filtering must return the same list when not asked to do anything.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ACG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter()
        self.assertEqual([read1, read2], list(result))

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

    def testFilterRemoveGaps(self):
        """
        Filtering must be able to remove gaps.
        """
        reads = Reads()
        reads.add(Read('id', '-AT--CG-'))
        result = reads.filter(removeGaps=True)
        self.assertEqual([Read('id', 'ATCG')], list(result))

    def testFilterNegativeRegex(self):
        """
        Filtering must be able to filter reads based on a negative regular
        expression.
        """
        reads = Reads()
        reads.add(Read('cats', 'ATCG'))
        reads.add(Read('kittens', 'ATCG'))
        reads.add(Read('dogs', 'ATCG'))
        reads.add(Read('puppies', 'ATCG'))
        reads.add(Read('lion', 'ATCG'))
        result = reads.filter(negativeTitleRegex='s')
        self.assertEqual([Read('lion', 'ATCG')], list(result))

    def testFilterPositiveRegex(self):
        """
        Filtering must be able to filter reads based on a positive regular
        expression.
        """
        reads = Reads()
        reads.add(Read('cats', 'ATCG'))
        reads.add(Read('kittens', 'ATCG'))
        reads.add(Read('dogs', 'ATCG'))
        reads.add(Read('puppies', 'ATCG'))
        reads.add(Read('lion', 'ATCG'))
        result = reads.filter(titleRegex='tt')
        self.assertEqual([Read('kittens', 'ATCG')], list(result))

    def testFilterWhitelist(self):
        """
        Filtering must be able to filter reads based on a whitelist.
        """
        reads = Reads()
        reads.add(Read('cats', 'ATCG'))
        reads.add(Read('kittens', 'ATCG'))
        reads.add(Read('dogs', 'ATCG'))
        reads.add(Read('puppies', 'ATCG'))
        reads.add(Read('lion', 'ATCG'))
        result = reads.filter(negativeTitleRegex='.', whitelist=['lion'])
        self.assertEqual([Read('lion', 'ATCG')], list(result))

    def testFilterBlacklist(self):
        """
        Filtering must be able to filter reads based on a blacklist.
        """
        reads = Reads()
        reads.add(Read('cats', 'ATCG'))
        reads.add(Read('kittens', 'ATCG'))
        reads.add(Read('dogs', 'ATCG'))
        reads.add(Read('puppies', 'ATCG'))
        reads.add(Read('lion', 'ATCG'))
        result = reads.filter(titleRegex='.',
                              blacklist=['cats', 'kittens', 'dogs', 'puppies'])
        self.assertEqual([Read('lion', 'ATCG')], list(result))

    def testFilterTruncateTitles(self):
        """
        Filtering must be able to filter reads based on a blacklist.
        """
        reads = Reads()
        reads.add(Read('cat 400', 'AA'))
        reads.add(Read('cat 500', 'GG'))
        result = reads.filter(truncateTitlesAfter='cat')
        self.assertEqual([Read('cat 400', 'AA')], list(result))

    def testFilterIndices(self):
        """
        Filtering must be able to filter reads based on their indices.
        """
        reads = Reads()
        reads.add(Read('cow', 'AA'))
        reads.add(Read('dog', 'GG'))
        reads.add(Read('cat', 'TT'))
        result = reads.filter(indices=set([1]))
        self.assertEqual([Read('dog', 'GG')], list(result))

    def testFilterHeadZero(self):
        """
        Filtering must be able to filter just the first N of a set of reads,
        including when N=0.
        """
        reads = Reads()
        reads.add(Read('cow', 'AA'))
        reads.add(Read('dog', 'GG'))
        reads.add(Read('cat', 'TT'))
        result = reads.filter(head=0)
        self.assertEqual([], list(result))

    def testFilterHead(self):
        """
        Filtering must be able to filter just the first N of a set of reads.
        """
        reads = Reads()
        reads.add(Read('cow', 'AA'))
        reads.add(Read('dog', 'GG'))
        reads.add(Read('cat', 'TT'))
        result = reads.filter(head=2)
        self.assertEqual([Read('cow', 'AA'), Read('dog', 'GG')], list(result))

    def testFilterDuplicates(self):
        """
        Filtering on sequence duplicates must work correctly.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(removeDuplicates=True)
        self.assertEqual([read1], list(result))


class TestSummarizePosition(TestCase):
    """
    Tests for the reads.summarizePosition function.
    """
    def testFrequenciesNoReads(self):
        """
        Must return empty counts if no reads are present.
        """
        reads = Reads()
        result = reads.summarizePosition(2)
        self.assertEqual({}, result['countAtPosition'])

    def testNumberOfExclusionsNoReads(self):
        """
        The excluded count must be zero if no reads are present.
        """
        reads = Reads()
        result = reads.summarizePosition(2)
        self.assertEqual(0, result['excludedCount'])

    def testExcludeShortSequences(self):
        """
        Sequences that are too short should be ignored.
        """
        reads = Reads()
        reads.add(Read('id1', 'agtcagtcagtc'))
        reads.add(Read('id2', 'acctg'))
        reads.add(Read('id3', 'atg'))
        result = reads.summarizePosition(9)
        self.assertEqual(2, result['excludedCount'])

    def testIndexLargerThanSequenceLength(self):
        """
        Must not count residues in sequences that are too short.
        """
        reads = Reads()
        reads.add(Read('id1', 'aaaaaa'))
        reads.add(Read('id2', 'aaca'))
        reads.add(Read('id3', 'aat'))
        result = reads.summarizePosition(5)
        self.assertEqual({'a': 1}, result['countAtPosition'])

    def testCorrectFrequencies(self):
        """
        Must return the correct frequencies.
        """
        reads = Reads()
        reads.add(Read('id1', 'aaaaaa'))
        reads.add(Read('id2', 'aata'))
        reads.add(Read('id3', 'aataaaaaa'))
        result = reads.summarizePosition(2)
        self.assertEqual({'a': 1, 't': 2}, result['countAtPosition'])

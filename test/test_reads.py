import six
from six.moves import builtins
from six import StringIO
from unittest import TestCase
from random import seed
from os import stat

try:
    from unittest.mock import patch, call
except ImportError:
    from mock import patch, call
from .mocking import mockOpen

from dark.aa import (
    BASIC_POSITIVE, HYDROPHOBIC, HYDROPHILIC, NEGATIVE, NONE, POLAR, SMALL,
    TINY)
from dark.fasta import FastaReads
from dark.hsp import HSP
from dark.reads import (
    Read, TranslatedRead, Reads, ReadsInRAM, DNARead, RNARead, DNAKozakRead,
    AARead, AAReadORF, AAReadWithX, SSAARead, SSAAReadWithX,
    readClassNameToClass)


class TestRead(TestCase):
    """
    Test the Read class.
    """
    def testGetitemReturnsNewRead(self):
        """
        __getitem__ must return a new Read instance.
        """
        self.assertIs(Read, Read('id', 'ACGT')[0:3].__class__)

    def testGetitemId(self):
        """
        __getitem__ must return a new Read instance with the same read id.
        """
        self.assertEqual('id-1234', Read('id-1234', 'ACGT')[0:3].id)

    def testGetitemSequence(self):
        """
        __getitem__ must return a Read instance with the expected sequence.
        """
        self.assertEqual('CG', Read('id', 'ACGT')[1:3].sequence)

    def testGetitemQuality(self):
        """
        __getitem__ must return a Read instance with the expected quality
        string.
        """
        self.assertEqual('12', Read('id', 'ACGT', '1234')[0:2].quality)

    def testGetitemLength(self):
        """
        __getitem__ must return a Read instance of the expected length.
        """
        self.assertEqual(3, len(Read('id-1234', 'ACGT')[0:3]))

    def testGetitemSingleIndex(self):
        """
        A single-index __getitem__ must return a length-one Read.
        """
        self.assertEqual(1, len(Read('id', 'ACGT')[0]))

    def testGetitemFullCopy(self):
        """
        A full copy __getitem__ must return the expected result.
        """
        self.assertEqual(Read('id', 'ACGT'),
                         Read('id', 'ACGT')[:])

    def testGetitemWithStep(self):
        """
        A stepped __getitem__ must return the expected result.
        """
        self.assertEqual(Read('id', 'AG', '13'),
                         Read('id', 'ACGT', '1234')[::2])

    def testGetitemReversed(self):
        """
        A reverse copy __getitem__ must return the expected result.
        """
        self.assertEqual(Read('id', 'TGCA', '4321'),
                         Read('id', 'ACGT', '1234')[::-1])

    def testUnequalLengths(self):
        """
        Attempting to construct a read whose sequence and quality strings are
        of different lengths must raise a ValueError.
        """
        error = ('Invalid read: sequence length \\(4\\) != quality '
                 'length \\(3\\)')
        with six.assertRaisesRegex(self, ValueError, error):
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
        error = "Format must be either 'fasta', 'fastq' or 'fasta-ss'\\."
        six.assertRaisesRegex(self, ValueError, error, read.toString,
                              'unknown')

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
        six.assertRaisesRegex(self, ValueError, error, read.toString, 'fastq')

    def testToFASTQ(self):
        """
        toString must return correct FASTA.
        """
        read = Read('id', 'ACGT', '!@#$')
        self.assertEqual('@id\nACGT\n+id\n!@#$\n', read.toString('fastq'))

    def testToDict(self):
        """
        toDict must return the correct dictionary.
        """
        read = Read('id3', 'ACGT', '!!2&')
        self.assertEqual(
            {
                'id': 'id3',
                'sequence': 'ACGT',
                'quality': '!!2&',
            },
            read.toDict())

    def testToDictNoQuality(self):
        """
        toDict must return the correct dictionary when the read has no quality
        string.
        """
        read = Read('id3', 'ACGT')
        self.assertEqual(
            {
                'id': 'id3',
                'sequence': 'ACGT',
                'quality': None,
            },
            read.toDict())

    def testFromDict(self):
        """
        fromDict must return the expected instance.
        """
        self.assertEqual(
            Read('id3', 'ACGT', '!!2&'),
            Read.fromDict({
                'id': 'id3',
                'sequence': 'ACGT',
                'quality': '!!2&',
            }))

    def testFromDictNoQuality(self):
        """
        fromDict must return the expected instance when the dictionary has no
        quality key.
        """
        self.assertEqual(
            Read('id3', 'ACGT'),
            Read.fromDict({
                'id': 'id3',
                'sequence': 'ACGT',
            }))

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

    def testHashViaSet(self):
        """
        If two identical reads are put into a set, the set must have size one.
        """
        read = Read('id', 'AA', '!!')
        self.assertEqual(1, len(set([read, read])))

    def testHashViaDict(self):
        """
        If two identical reads are used as keys in a dict, the dict must have
        size one.
        """
        read = Read('id', 'AA', '!!')
        self.assertEqual(1, len(dict.fromkeys([read, read])))

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

    def testWalkHSPLeftOverhangingMatchNoWhiskers(self):
        """
        If the HSP specifies that the entire read matches the subject, and
        also extends to the left of the subject, walkHSP must return the
        correct results when it is told to not include whiskers.

        Subject:       GT.....
        Read:        ACGT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=2, readEnd=4, readStartInSubject=-2,
                  readEndInSubject=2, subjectStart=0, subjectEnd=2,
                  readMatchedSequence='GT', subjectMatchedSequence='GT')
        self.assertEqual([(0, 'G', True),
                          (1, 'T', True)],
                         list(read.walkHSP(hsp, includeWhiskers=False)))

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

    def testWalkHSPRightOverhangingMatchNoWhiskers(self):
        """
        If the HSP specifies that the entire read matches the subject, and
        also extends to the right of the subject, walkHSP must return the
        correct results when it is told to not include whiskers.

        Subject:       AC
        Read:          ACGT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=0, readEnd=2, readStartInSubject=10,
                  readEndInSubject=14, subjectStart=10, subjectEnd=12,
                  readMatchedSequence='AC', subjectMatchedSequence='AC')
        self.assertEqual([(10, 'A', True),
                          (11, 'C', True)],
                         list(read.walkHSP(hsp, includeWhiskers=False)))

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

    def testWalkHSPLeftAndRightOverhangingMatchNoWhiskers(self):
        """
        If the HSP specifies that the read matches the entire subject, and
        also extends to both the left and right of the subject, walkHSP must
        return the correct results when it is told to not include whiskers.

        Subject:        CG
        Read:          ACGT
        """
        read = Read('id', 'ACGT')
        hsp = HSP(33, readStart=1, readEnd=3, readStartInSubject=10,
                  readEndInSubject=14, subjectStart=11, subjectEnd=13,
                  readMatchedSequence='CG', subjectMatchedSequence='CG')
        self.assertEqual([(11, 'C', True),
                          (12, 'G', True)],
                         list(read.walkHSP(hsp, includeWhiskers=False)))

    def testCheckAlphabetwithReadMustBePermissive(self):
        """
        The checkAlphabet function must return the expected alphabet if a
        dark.Read is passed.
        """
        read = Read('id', 'ARSTGATGCASASASASASAS')
        self.assertEqual(set('ACGSRT'), read.checkAlphabet())

    def testCheckAlphabetAAReadMatchingReturnTrue(self):
        """
        If an AA read with an AARead readClass is passed in, the checkAlphabet
        function must return the alphabet of the sequence.
        """
        read = AARead('id', 'ARSTGATGCASASASASASAS')
        self.assertEqual(set('ACGSRT'), read.checkAlphabet())

    def testCheckAlphabetDNAReadMatchingReturnTrue(self):
        """
        If a DNA read with a DNARead readClass is passed in, the checkAlphabet
        function must return the alphabet of the sequence.
        """
        read = DNARead('id', 'AAATTAACGGGCCTAGG')
        self.assertEqual(set('ACTG'), read.checkAlphabet())

    def testCheckAlphabetAAReadNotMatchingRaise(self):
        """
        If an AA read with a DNARead readClass is passed in, the checkAlphabet
        function must raise an IndexError.
        """
        read = AARead('id', 'AAATTAACGGGCCTAGG')
        error = "It looks like a DNA sequence has been passed to AARead()."
        six.assertRaisesRegex(self, ValueError, error, read.checkAlphabet)

    def testCheckAlphabetDNAReadNotMatchingRaise(self):
        """
        If a DNA read with an AARead readClass is passed in, the checkAlphabet
        function must raise an IndexError.
        """
        read = DNARead('id', 'ARSTGATGCASASASASASAS')
        error = ("Read alphabet \\('ACGRST'\\) is not a subset of expected "
                 "alphabet \\('ACGT'\\) for read class DNARead.")
        six.assertRaisesRegex(self, ValueError, error, read.checkAlphabet)

    def testKeepSites(self):
        """
        If only a certain set of sites should be kept, newFromSites should
        return a read with the correct sequence.
        """
        self.assertEqual(Read('id1', 'TC'),
                         Read('id1', 'ATCGAT').newFromSites({1, 2}))

    def testKeepSitesNoSites(self):
        """
        If only the empty set of sites should be kept, newFromSites should
        return a read with the correct (empty) sequence.
        """
        self.assertEqual(Read('id1', ''),
                         Read('id1', 'ATCGAT').newFromSites(set()))

    def testKeepSitesAllSites(self):
        """
        If all sites should be kept, newFromSites should return a read with
        the correct (full) sequence.
        """
        self.assertEqual(Read('id1', 'ATCGAT'),
                         Read('id1', 'ATCGAT').newFromSites(set(range(6))))

    def testKeepSitesWithQuality(self):
        """
        If only a certain set of sites should be kept, newFromSites should
        return a read with the correct sequence and quality.
        """
        self.assertEqual(Read('id1', 'TC', '23'),
                         Read('id1', 'ATCGAT', '123456').newFromSites(
                             {1, 2}))

    def testKeepSitesOutOfRange(self):
        """
        If only a certain set of sites should be kept, but the kept sites
        are higher than the length of the input sequences, newFromSites
        should return a read with the correct (empty) sequence.
        """
        self.assertEqual(Read('id1', ''),
                         Read('id1', 'ATCGAT').newFromSites({100, 200}))

    def testRemoveSites(self):
        """
        If only a certain set of sites should be removed, newFromSites
        should return a read with the correct sequence.
        """
        self.assertEqual(Read('id1', 'AGAT'),
                         Read('id1', 'ATCGAT').newFromSites({1, 2},
                                                            exclude=True))

    def testRemoveSitesNoSites(self):
        """
        If no sites should be removed, newFromSites should return a read
        with the correct (full) sequence.
        """
        self.assertEqual(Read('id1', 'ATCGAT'),
                         Read('id1', 'ATCGAT').newFromSites(set(),
                                                            exclude=True))

    def testRemoveSitesAllSites(self):
        """
        If all sites should be removed, newFromSites should return a read
        with the correct (empty) sequence.
        """
        self.assertEqual(Read('id1', ''),
                         Read('id1', 'ATCGAT').newFromSites(set(range(6)),
                                                            exclude=True))

    def testRemoveSitesWithQuality(self):
        """
        If only a certain set of sites should be removed, newFromSites
        should return a read with the correct sequence and quality.
        """
        self.assertEqual(Read('id1', 'AGAT', '1456'),
                         Read('id1', 'ATCGAT', '123456').newFromSites(
                             {1, 2}, exclude=True))

    def testRemoveSitesOutOfRange(self):
        """
        If only a certain set of sites should be removed, but the removed
        sites are higher than the length of the input sequences,
        newFromSites should return a read with the correct (full) sequence.
        """
        self.assertEqual(Read('id1', 'ATCGAT'),
                         Read('id1', 'ATCGAT').newFromSites({100, 200},
                                                            exclude=True))

    def testReverseNoQuality(self):
        """
        The reverse method must work as expected when there is no quality
        string.
        """
        self.assertEqual(Read('id1', 'ATCGAT'),
                         Read('id1', 'TAGCTA').reverse())

    def testReverseWithQuality(self):
        """
        The reverse method must work as expected when there is a quality
        string.
        """
        self.assertEqual(Read('id1', 'ATCGAT', '123456'),
                         Read('id1', 'TAGCTA', '654321').reverse())


class TestDNARead(TestCase):
    """
    Tests for the DNARead class.
    """
    def testGetitemReturnsNewDNARead(self):
        """
        __getitem__ must return a new DNARead instance.
        """
        self.assertIs(DNARead, DNARead('id', 'ACGT')[0:3].__class__)

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

    def testReverseComplementLowercaseLetters(self):
        """
        The reverseComplement function must correctly reverse complement
        lowercase letters. The issue is described here:
        https://github.com/acorg/dark-matter/issues/662
        """
        read = DNARead('id', 'CAGCAGctgcagcaccagcaccagcagcttcCACAT')
        expected = ('ATGTGgaagctgctggtgctggtgctgcagCTGCTG')
        self.assertEqual(expected, read.reverseComplement().sequence)

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
    def testGetitemReturnsNewRNARead(self):
        """
        __getitem__ must return a new RNARead instance.
        """
        self.assertIs(RNARead, RNARead('id', 'ACGU')[0:3].__class__)

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


class TestDNAKozakRead(TestCase):
    """
    Test the DNAKozakRead class.
    """
    def testSequence(self):
        """
        A DNAKozakRead instance must have the correct sequence.
        """
        originalRead = DNARead('id', 'AAGTAAGGGCTGTGA')
        read = DNAKozakRead(originalRead, 0, 4, 100.0)
        self.assertEqual('AAGT', read.sequence)

    def testStart(self):
        """
        A DNAKozakRead instance must store the correct start offset.
        """
        originalRead = DNARead('id', 'AAGTAAGGGCTGTGA')
        read = DNAKozakRead(originalRead, 2, 4, 100.0)
        self.assertEqual(2, read.start)

    def testStop(self):
        """
        A DNAKozakRead instance must store the correct stop offset.
        """
        originalRead = DNARead('id', 'AAGTAAGGGCTGTGA')
        read = DNAKozakRead(originalRead, 2, 4, 100.0)
        self.assertEqual(4, read.stop)
        read = DNAKozakRead(originalRead, 4, 10, 100.0)
        self.assertEqual(10, read.stop)

    def testStartNegative(self):
        """
        A DNAKozakRead start offset must not be less than zero.
        """
        originalRead = DNARead('id', 'AAGTAAGGGCTGTGA')
        error = r'^start offset \(-1\) less than zero$'
        six.assertRaisesRegex(
            self, ValueError, error, DNAKozakRead, originalRead, -1, 6, 100.0)

    def testStartGreaterThanStop(self):
        """
        A DNAKozakRead start offset must not be greater than its stop offset.
        """
        originalRead = DNARead('id', 'AAGTAAGGGCTGTGA')
        error = 'start offset \\(4\\) greater than stop offset \\(0\\)'
        six.assertRaisesRegex(
            self, ValueError, error, DNAKozakRead, originalRead, 4, 0, 100.0)

    def testStopGreaterThanOriginalSequenceLength(self):
        """
        A DNAKozakRead stop offset must not be greater than the length of the
        original sequence.
        """
        originalRead = DNARead('id', 'AAGTAA')
        error = 'stop offset \\(10\\) > original read length \\(6\\)'
        six.assertRaisesRegex(
            self, ValueError, error, DNAKozakRead, originalRead, 0, 10, 100.0)

    def testEqualFunction(self):
        """
        A DNAKozakRead needs to compare correctly to an equal other
        DNAKozakRead.
        """
        originalRead = DNARead('id', 'AAGTAAGGGCTGTGA')
        kozakRead1 = DNAKozakRead(originalRead, 0, 4, 100.0)
        kozakRead2 = DNAKozakRead(originalRead, 0, 4, 100.0)
        self.assertEqual(kozakRead1, kozakRead2)

    def testEqualFunctionDifferentOriginalSequence(self):
        """
        A DNAKozakRead needs to compare correctly to a different other
        DNAKozakRead.
        """
        originalRead1 = DNARead('id', 'AAGTAAGGGCTGTGA')
        originalRead2 = DNARead('id', 'AAGTAAGGGCTGTGAAA')
        kozakRead1 = DNAKozakRead(originalRead1, 0, 4, 100.0)
        kozakRead2 = DNAKozakRead(originalRead2, 0, 4, 100.0)
        self.assertNotEqual(kozakRead1, kozakRead2)

    def testEqualFunctionDifferentKozakSequence(self):
        """
        A DNAKozakRead needs to compare correctly to a different other
        DNAKozakRead.
        """
        originalRead1 = DNARead('id', 'AAGTAAGGGCTGTGA')
        originalRead2 = DNARead('id', 'AAGTAAGGGCTGTGAAA')
        kozakRead1 = DNAKozakRead(originalRead1, 0, 4, 100.0)
        kozakRead2 = DNAKozakRead(originalRead2, 1, 5, 100.0)
        self.assertNotEqual(kozakRead1, kozakRead2)


class TestAARead(TestCase):
    """
    Tests for the AARead class.
    """
    def testGetitemReturnsNewAARead(self):
        """
        __getitem__ must return a new AARead instance.
        """
        self.assertIs(AARead, AARead('id', 'ACGU')[0:3].__class__)

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
                    'polar requirement': -0.463414634146,
                    'hydropathy': 0.4,
                    'iep': -0.191489361702,
                    'hydroxythiolation': -0.265160523187,
                    'aromaticity': -0.550128534704,
                    'hydrogenation': 0.8973042362,
                    'composition': -1.0
                },
                {
                    'polarity': 1.0,
                    'aliphaticity': -0.818181818182,
                    'volume': -0.389221556886,
                    'polar requirement': 1.0,
                    'hydropathy': -0.777777777778,
                    'iep': -1.0,
                    'hydroxythiolation': -0.348394768133,
                    'aromaticity': -1.0,
                    'hydrogenation': -0.90243902439,
                    'composition': 0.00363636363636
                },
                {
                    'polarity': 0.382716049383,
                    'aliphaticity': -0.157024793388,
                    'volume': 0.449101796407,
                    'polar requirement': 0.0487804878049,
                    'hydropathy': -1.0,
                    'iep': 1.0,
                    'hydroxythiolation': -0.51486325802,
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
        orfs = list(read.ORFs(True))
        self.assertEqual(0, len(orfs))

    def testORFsEmptySequenceWithStartStop(self):
        """
        An AA read with just a start and stop codon must not have any ORFs.
        """
        read = AARead('id', 'M*')
        orfs = list(read.ORFs(False))
        self.assertEqual(0, len(orfs))

    def testORFsEmptySequenceWithStartStopOpenORFs(self):
        """
        An AA read with just a start and stop codon must have one ORF.
        """
        read = AARead('id', 'M*')
        orfs = list(read.ORFs(True))
        self.assertEqual(1, len(orfs))

    def testORFsEmptySequenceWithStopStartOpenORFs(self):
        """
        An AA read with just a start and stop codon must have one ORF.
        """
        read = AARead('id', '*M')
        orfs = list(read.ORFs(True))
        self.assertEqual(0, len(orfs))

    def testORFsEmptySequenceWithStart(self):
        """
        An AA read with just a start codon must not have any ORFs.
        """
        read = AARead('id', 'M')
        orfs = list(read.ORFs(False))
        self.assertEqual(0, len(orfs))

    def testORFsEmptySequenceWithStartOpenORFs(self):
        """
        An AA read with just a start codon must have one ORF.
        """
        read = AARead('id', 'M')
        orfs = list(read.ORFs(True))
        self.assertEqual(1, len(orfs))

    def testORFsSequenceWithOneAAOpenORFs(self):
        """
        An AA read with just a start codon must have one ORF.
        """
        read = AARead('id', 'A')
        orfs = list(read.ORFs(True))
        self.assertEqual(1, len(orfs))

    def testORFsWithOneStopCodon(self):
        """
        An AA read of a single stop codon must not have any ORFs.
        """
        read = AARead('id', '*')
        orfs = list(read.ORFs(False))
        self.assertEqual(0, len(orfs))

    def testORFsWithOneStopCodonOpenORFs(self):
        """
        An AA read of a single stop codon must not have any ORFs.
        """
        read = AARead('id', '*')
        orfs = list(read.ORFs(True))
        self.assertEqual(0, len(orfs))

    def testORFsWithTwoStopCodons(self):
        """
        An AA read of two stop codons must not have any ORFs.
        """
        read = AARead('id', '**')
        orfs = list(read.ORFs(True))
        self.assertEqual(0, len(orfs))

    def testORFsWithJustStartsAndStops(self):
        """
        An AA read of only start and stop codons must not have any ORFs.
        """
        read = AARead('id', '**MM*M**MMM*')
        orfs = list(read.ORFs(False))
        self.assertEqual(2, len(orfs))

    def testOpenOpenORF(self):
        """
        An AA read that contains no start or stop codons should result in
        just one AAReadORF when its ORFs method is called. The ORF must have
        the correct start/stop offsets and its left and right side must be
        marked as open.
        """
        read = AARead('id', 'ADRADR')
        orfs = list(read.ORFs(True))
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
        orfs = list(read.ORFs(True))
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
        orfs = list(read.ORFs(True))
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(6, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:6]')

    def testCloseOpenORFWithMultipleStarts(self):
        """
        An AA read that contains multiple initial start codons but no
        stop codon should result in just one AAReadORF when its ORFs method
        is called. The ORF must have the correct start/stop offsets and its
        left and right sides must be marked as closed and open,
        respectively.
        """
        read = AARead('id', 'MMMADRADR')
        orfs = list(read.ORFs(True))
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('MMMADRADR', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(9, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:9)')

    def testCloseCloseORF(self):
        """
        An AA read that contains a start and a stop codon should result in
        just one AAReadORF when its ORFs method is called. The ORF must have
        the correct start/stop offsets and its left and right sides must be
        both marked as closed.
        """
        read = AARead('id', 'MADRADR*')
        orfs = list(read.ORFs(False))
        self.assertEqual(1, len(orfs))
        orf = orfs[0]
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(1, orf.start)
        self.assertEqual(7, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[1:7]')

    def testCloseCloseORFWithJunk(self):
        """
        An AA read that contains a start and a stop codon should result in
        just one AAReadORF when its ORFs method is called. The ORF must have
        the correct start/stop offsets and its left and right sides must be
        both marked as closed.
        """
        read = AARead('id', 'AAAMADRADR*')
        [orf] = list(read.ORFs(False))
        self.assertEqual('ADRADR', orf.sequence)
        self.assertEqual(4, orf.start)
        self.assertEqual(10, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[4:10]')

    def testOpenCloseThenCloseOpenORF(self):
        """
        An AA read that contains an ORF that is left open and right closed
        followed by an ORF that is left closed and right open must have the
        ORFs detected correctly when its ORFs method is called.
        """
        read = AARead('id', 'ADR*MRRR')
        orfs = list(read.ORFs(True))
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
        read = AARead('id', '*MADR*MRRR')
        orfs = list(read.ORFs(True))
        self.assertEqual(2, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(2, orf.start)
        self.assertEqual(5, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[2:5]')

        orf = orfs[1]
        self.assertEqual('RRR', orf.sequence)
        self.assertEqual(7, orf.start)
        self.assertEqual(10, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertTrue(orf.openRight)
        self.assertEqual(orf.id, 'id-[7:10)')

    def testCloseCloseThenCloseCloseORF(self):
        """
        An AA read that contains two ORFs that are both left and right closed
        must have the ORFs detected correctly when its ORFs method is called.
        """
        read = AARead('id', 'MADR*MRRR*')
        orfs = list(read.ORFs(False))
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
        orfs = list(read.ORFs(True))
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

    def testCloseCloseThenCloseCloseThenNothingORF(self):
        """
        An AA read that contains an ORF that is left and right closed
        followed by an internal ORF, followed by an ORF that is left closed
        and right open must have the ORFs detected correctly when its ORFs
        method is called.
        """
        read = AARead('id', 'MADR*MAAA*MRRR')
        orfs = list(read.ORFs(False))
        self.assertEqual(2, len(orfs))

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

    def testCloseCloseThenCloseCloseThenCloseCloseORF(self):
        """
        An AA read that contains an ORF that is left and right closed
        followed by an internal ORF, followed by an ORF that is left and
        right closed must have the ORFs detected correctly when its ORFs
        method is called.
        """
        read = AARead('id', 'MADR*MAAA*MRRR*')
        orfs = list(read.ORFs(False))
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
        orfs = list(read.ORFs(True))
        self.assertEqual(4, len(orfs))

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
        self.assertEqual('MM', orf.sequence)
        self.assertEqual(16, orf.start)
        self.assertEqual(18, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[16:18]')

        orf = orfs[3]
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
        orfs = list(read.ORFs(True))
        self.assertEqual(4, len(orfs))

        orf = orfs[0]
        self.assertEqual('ADR', orf.sequence)
        self.assertEqual(3, orf.start)
        self.assertEqual(6, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[3:6]')

        orf = orfs[1]
        self.assertEqual('M', orf.sequence)
        self.assertEqual(10, orf.start)
        self.assertEqual(11, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[10:11]')

        orf = orfs[2]
        self.assertEqual('AAA', orf.sequence)
        self.assertEqual(14, orf.start)
        self.assertEqual(17, orf.stop)
        self.assertFalse(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-[14:17]')

        orf = orfs[3]
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
        orfs = list(read.ORFs(False))
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
        orfs = list(read.ORFs(True))
        self.assertEqual(1, len(orfs))

        orf = orfs[0]
        self.assertEqual('KK', orf.sequence)
        self.assertEqual(0, orf.start)
        self.assertEqual(2, orf.stop)
        self.assertTrue(orf.openLeft)
        self.assertFalse(orf.openRight)
        self.assertEqual(orf.id, 'id-(0:2]')


class TestAAReadWithX(TestCase):
    """
    Tests for the TestAAReadWithX class.
    """
    def testGetitemReturnsNewAAReadWithX(self):
        """
        __getitem__ must return a new AAReadWithX instance.
        """
        self.assertIs(AAReadWithX, AAReadWithX('id', 'ACGU')[0:3].__class__)

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
        error = 'start offset \\(4\\) greater than stop offset \\(0\\)'
        six.assertRaisesRegex(
            self, ValueError, error, AAReadORF, originalRead, 4, 0, True, True)

    def testStartNegative(self):
        """
        An AAReadORF start offset must not be less than zero.
        """
        originalRead = AARead('id', 'ADRADR')
        error = 'start offset \\(-1\\) less than zero'
        six.assertRaisesRegex(
            self, ValueError, error, AAReadORF, originalRead, -1, 6, True,
            True)

    def testStopGreaterThanOriginalSequenceLength(self):
        """
        An AAReadORF stop offset must not be greater than the length of the
        original sequence.
        """
        originalRead = AARead('id', 'ADRADR')
        error = 'stop offset \\(10\\) > original read length \\(6\\)'
        six.assertRaisesRegex(
            self, ValueError, error, AAReadORF, originalRead, 0, 10, True,
            True)

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

    def testToDict(self):
        """
        toDict must return the correct dictionary.
        """
        originalRead = AARead('id3', 'ACGT', '!!2&')
        read = AAReadORF(originalRead, 1, 3, True, False)
        self.assertEqual(
            {
                'id': 'id3-(1:3]',
                'sequence': 'CG',
                'quality': '!2',
                'start': 1,
                'stop': 3,
                'openLeft': True,
                'openRight': False,
            },
            read.toDict())

    def testFromDict(self):
        """
        fromDict must return the expected instance.
        """
        originalRead = AARead('id3', 'ACGT', '!!2&')
        self.assertEqual(
            AAReadORF(originalRead, 1, 3, True, False),
            AAReadORF.fromDict({
                'id': 'id3-(1:3]',
                'sequence': 'CG',
                'quality': '!2',
                'start': 1,
                'stop': 3,
                'openLeft': True,
                'openRight': False,
            }))


class _TestSSAAReadMixin(object):
    """
    Mixin class with tests for the SSAARead and SSAAReadWithX classes.
    """
    def testSequenceLengthMatchesStructureLength(self):
        """
        An SSAARead (or SSAAReadWithX) must have sequence and structure
        lengths that are the same.
        """
        error = (
            '^Invalid read: sequence length \\(4\\) != structure '
            'length \\(3\\)$')
        with six.assertRaisesRegex(self, ValueError, error):
            self.CLASS('id', 'ACGT', '!!!')

    def testCorrectAttributes(self):
        """
        An SSAARead (or SSAAReadWithX) must have the correct attributes.
        """
        read = self.CLASS('id', 'AFGGCED', 'HHH  HH')
        self.assertEqual('id', read.id)
        self.assertEqual('AFGGCED', read.sequence)
        self.assertEqual('HHH  HH', read.structure)
        self.assertIs(None, read.quality)

    def testReads(self):
        """
        It must be possible to make a dark.Reads object out of SSAARead
        (or SSAAReadWithX) instances, and the result must have the correct
        length.
        """
        reads = Reads()
        reads.add(self.CLASS('id1', 'AFGGCED', 'HHHHHHH'))
        reads.add(self.CLASS('id2', 'AFGGKLL', 'HHHHIII'))
        self.assertEqual(2, len(list(reads)))

    def testGetitemReturnsNewRead(self):
        """
        __getitem__ must return a new instance of the correct class.
        """
        self.assertIs(self.CLASS,
                      self.CLASS('id', 'ACGT', 'HHHH')[0:3].__class__)

    def testGetitemId(self):
        """
        __getitem__ must return a new instance with the same read id.
        """
        self.assertEqual('id-12', self.CLASS('id-12', 'FFRR', 'HHHH')[0:3].id)

    def testGetitemSequence(self):
        """
        __getitem__ must return a new instance with the expected sequence.
        """
        self.assertEqual('RM', self.CLASS('id', 'FRML', 'HHHH')[1:3].sequence)

    def testGetitemStructure(self):
        """
        __getitem__ must return a new instance with the expected structure.
        """
        self.assertEqual('HE', self.CLASS('id', 'FFRR', 'HEIS')[0:2].structure)

    def testGetitemLength(self):
        """
        __getitem__ must return a new instance of the expected length.
        """
        self.assertEqual(3, len(self.CLASS('id-1234', 'FFMM', 'HHHH')[0:3]))

    def testGetitemSingleIndex(self):
        """
        A single-index __getitem__ must return a length-one instance.
        """
        self.assertEqual(1, len(self.CLASS('id', 'FRFR', 'HHHH')[0]))

    def testGetitemFullCopy(self):
        """
        A full copy __getitem__ must return the expected result.
        """
        self.assertEqual(self.CLASS('id', 'RRFF', 'HEIH'),
                         self.CLASS('id', 'RRFF', 'HEIH')[:])

    def testGetitemWithStep(self):
        """
        A stepped __getitem__ must return the expected result.
        """
        self.assertEqual(self.CLASS('id', 'FR', 'HS'),
                         self.CLASS('id', 'FMRL', 'HESE')[::2])

    def testGetitemReversed(self):
        """
        A reverse copy __getitem__ must return the expected result.
        """
        self.assertEqual(self.CLASS('id', 'FRML', 'HESB'),
                         self.CLASS('id', 'LMRF', 'BSEH')[::-1])

    def testToString(self):
        """
        toString must return the expected 2 FASTA records (one of which is
        the structure information).
        """
        self.assertEqual(
            '>id-1234\n'
            'FFMM\n'
            '>id-1234:structure\n'
            'HHHH\n',
            self.CLASS('id-1234', 'FFMM', 'HHHH').toString())

    def testToStringWithStructureSuffix(self):
        """
        toString must return the expected 2 FASTA records when given a
        specific structure id suffix.
        """
        self.assertEqual(
            '>id-12\n'
            'FFMM\n'
            '>id-12:x\n'
            'HHHH\n',
            self.CLASS('id-12', 'FFMM', 'HHHH').toString(structureSuffix=':x'))

    def testToStringWithExplicitFastaSSFormat(self):
        """
        toString must return the expected 2 FASTA records when 'fasta-ss' is
        passed as the C{format_} argument.
        """
        self.assertEqual(
            '>id-1234\n'
            'FFMM\n'
            '>id-1234:structure\n'
            'HHHH\n',
            self.CLASS('id-1234', 'FFMM', 'HHHH').toString(format_='fasta-ss'))

    def testToStringWithExplicitFastaFormat(self):
        """
        toString must return normal FASTA when 'fasta' is passed as the
        C{format_} argument.
        """
        self.assertEqual(
            '>id-1234\n'
            'FFMM\n',
            self.CLASS('id-1234', 'FFMM', 'HHHH').toString(format_='fasta'))

    def testToStringWithUnknownFormat(self):
        """
        toString must raise ValueError when something other than 'fasta' or
        'fasta-ss' is passed as the C{format_} argument.
        """
        read = self.CLASS('id-1234', 'FFMM', 'HHHH')
        error = "^Format must be either 'fasta', 'fastq' or 'fasta-ss'\\."
        six.assertRaisesRegex(
            self, ValueError, error, read.toString, format_='pasta')

    def testToDict(self):
        """
        toDict must return the correct dictionary.
        """
        read = SSAARead('id3', 'ACGT', 'HHEE')
        self.assertEqual(
            {
                'id': 'id3',
                'sequence': 'ACGT',
                'structure': 'HHEE',
            },
            read.toDict())

    def testFromDict(self):
        """
        fromDict must return the expected instance.
        """
        self.assertEqual(
            SSAARead('id3', 'ACGT', 'HHEE'),
            SSAARead.fromDict({
                'id': 'id3',
                'sequence': 'ACGT',
                'structure': 'HHEE',
            }))

    def testHashDiffersIfIdDiffers(self):
        """
        The __hash__ value for two reads must differ if their ids differ.
        """
        self.assertNotEqual(hash(self.CLASS('id1', 'AA', 'HH')),
                            hash(self.CLASS('id2', 'AA', 'HH')))

    def testHashDiffersIfSequenceDiffers(self):
        """
        The __hash__ value for two reads must differ if their sequence strings
        differ.
        """
        self.assertNotEqual(hash(self.CLASS('id', 'MMR', 'HHH')),
                            hash(self.CLASS('id', 'MMF', 'HHH')))

    def testHashDiffersIfStructureDiffers(self):
        """
        The __hash__ value for two reads must differ if their structure strings
        differ.
        """
        self.assertNotEqual(hash(self.CLASS('id', 'AA', 'HH')),
                            hash(self.CLASS('id', 'AA', 'HE')))

    def testHashViaSet(self):
        """
        If two identical reads are put into a set, the set must have size one.
        """
        read = self.CLASS('id', 'AA', 'HH')
        self.assertEqual(1, len(set([read, read])))

    def testHashViaDict(self):
        """
        If two identical reads are used as keys in a dict, the dict must have
        size one.
        """
        read = self.CLASS('id', 'AA', 'HH')
        self.assertEqual(1, len(dict.fromkeys([read, read])))

    def testKeepSites(self):
        """
        If only a certain set of sites should be kept, newFromSites should
        return a read with the correct sequence.
        """
        self.assertEqual(self.CLASS('id1', 'TC', 'HG'),
                         self.CLASS('id1', 'ATCGAT', 'CHGCCX').newFromSites(
                             {1, 2}))

    def testKeepSitesNoSites(self):
        """
        If an empty set of sites should be kept, newFromSites should
        return a read with the correct (empty) sequence.
        """
        self.assertEqual(self.CLASS('id1', '', ''),
                         self.CLASS('id1', 'ATCGAT', 'CHGCCC').newFromSites(
                             set()))

    def testKeepSitesAllSites(self):
        """
        If all sites should be kept, newFromSites should return a read
        with the correct (full) sequence.
        """
        self.assertEqual(self.CLASS('id1', 'ATCGAT', 'CHGCCC'),
                         self.CLASS('id1', 'ATCGAT', 'CHGCCC').newFromSites(
                             set(range(6))))

    def testKeepSitesOutOfRange(self):
        """
        If only a certain set of sites should be kept, but the kept sites
        are higher than the length of the input sequences, newFromSites
        should return a read with the correct (empty) sequence.
        """
        self.assertEqual(self.CLASS('id1', '', ''),
                         self.CLASS('id1', 'ATCGAT', 'CHGCCC').newFromSites(
                             {100, 200}))

    def testRemoveSites(self):
        """
        If only a certain set of sites should be removed, newFromSites
        should return a read with the correct sequence.
        """
        self.assertEqual(self.CLASS('id1', 'AGAT', 'CSHG'),
                         self.CLASS('id1', 'ATCGAT', 'CXXSHG').newFromSites(
                             {1, 2}, exclude=True))

    def testRemoveSitesNoSites(self):
        """
        If the empty set of sites should be removed, newFromSites
        should return a read with the correct (full) sequence.
        """
        self.assertEqual(self.CLASS('id1', 'ATCGAT', 'CXXSHG'),
                         self.CLASS('id1', 'ATCGAT', 'CXXSHG').newFromSites(
                             set(), exclude=True))

    def testRemoveSitesAllSites(self):
        """
        If all sites should be removed, newFromSites should return a read
        with the correct (empty) sequence.
        """
        self.assertEqual(self.CLASS('id1', '', ''),
                         self.CLASS('id1', 'ATCGAT', 'CXXSHG').newFromSites(
                             set(range(6)), exclude=True))

    def testRemoveSitesOutOfRange(self):
        """
        If only a certain set of sites should be removed, but the removed
        sites are higher than the length of the input sequences,
        newFromSites should return a read with the correct (full) sequence.
        """
        self.assertEqual(self.CLASS('id1', 'ATCGAT', 'HGSTCC'),
                         self.CLASS('id1', 'ATCGAT', 'HGSTCC').newFromSites(
                             {100, 200}, exclude=True))


class TestSSAARead(TestCase, _TestSSAAReadMixin):
    """
    Tests for the SSAARead class.
    """
    CLASS = SSAARead


class TestSSAAReadWithX(TestCase, _TestSSAAReadMixin):
    """
    Tests for the SSAAReadWithX class.
    """
    CLASS = SSAAReadWithX

    def testSequenceContainingX(self):
        """
        An SSAAReadWithX must be able to contain an 'X' character.
        """
        self.assertEqual('AFGX', SSAAReadWithX('id', 'AFGX', 'HHHH').sequence)


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
        six.assertRaisesRegex(self, ValueError, error, TranslatedRead, read,
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
        read = Read('id', 'acctagatggttgtttag')
        translated = TranslatedRead(read, 'T*MVV*', 0)
        self.assertEqual(2, translated.maximumORFLength())

    def testMaximumORFLengthNoOpenORF(self):
        """
        The maximumORFLength function must return the correct value if
        open ORFs are not allowed.
        """
        read = Read('id', 'atgacctagatggttgtttag')
        translated = TranslatedRead(read, 'MT*MVV*', 0)
        self.assertEqual(2, translated.maximumORFLength(False))

    def testToDict(self):
        """
        toDict must return the correct dictionary.
        """
        originalRead = AARead('id3', 'ACGT', '!!2&')
        read = TranslatedRead(originalRead, 'MMMM', 0, True)
        self.assertEqual(
            {
                'id': 'id3-frame0rc',
                'sequence': 'MMMM',
                'quality': None,
                'frame': 0,
                'reverseComplemented': True,
            },
            read.toDict())

    def testFromDict(self):
        """
        fromDict must return the expected instance.
        """
        originalRead = AARead('id3', 'ACGT')
        self.assertEqual(
            TranslatedRead(originalRead, 'MMMM', 0, True),
            TranslatedRead.fromDict({
                'id': 'id3-frame0rc',
                'sequence': 'MMMM',
                'quality': None,
                'frame': 0,
                'reverseComplemented': True,
            }))


class TestReadClassNameToClass(TestCase):
    """
    Test that the light.reads.readClassNameToClass dictionary is correct.
    """
    def testNames(self):
        self.assertEqual(9, len(readClassNameToClass))
        self.assertIs(AARead, readClassNameToClass['AARead'])
        self.assertIs(AAReadORF, readClassNameToClass['AAReadORF'])
        self.assertIs(AAReadWithX, readClassNameToClass['AAReadWithX'])
        self.assertIs(DNARead, readClassNameToClass['DNARead'])
        self.assertIs(RNARead, readClassNameToClass['RNARead'])
        self.assertIs(Read, readClassNameToClass['Read'])
        self.assertIs(SSAARead, readClassNameToClass['SSAARead'])
        self.assertIs(SSAAReadWithX, readClassNameToClass['SSAAReadWithX'])
        self.assertIs(TranslatedRead, readClassNameToClass['TranslatedRead'])


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
        self.assertEqual(0, len(list(reads)))

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

    def testEmptyInitialReads(self):
        """
        A Reads instance must be able to accept an empty initial iterable of
        reads.
        """
        reads = Reads([])
        self.assertEqual([], list(reads))

    def testInitialReads(self):
        """
        A Reads instance must be able to accept a non-empty initial iterable
        of reads.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads = Reads(initialReads=[read1, read2])
        self.assertEqual([read1, read2], list(reads))

    def testManuallyAddedReadsLength(self):
        """
        A Reads instance with reads added manually must have the correct
        length.
        """
        reads = Reads()
        reads.add(Read('id1', 'AT'))
        reads.add(Read('id2', 'AC'))
        self.assertEqual(2, len(list(reads)))

    def testSubclass(self):
        """
        A Reads subclass with an iter method must result in an instance
        with a correct iterator.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')

        class ReadsSubclass(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = ReadsSubclass()
        self.assertEqual([read1, read2], list(reads))

    def testSubclassLength(self):
        """
        A Reads subclass with an iter method must result in an instance
        with a correct length.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')

        class ReadsSubclass(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = ReadsSubclass()
        self.assertEqual(2, len(list(reads)))

    def testRepeatedIter(self):
        """
        A Reads subclass with an iter method must be able to be listed
        more than once.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')

        class ReadsSubclass(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = ReadsSubclass()
        self.assertEqual([read1, read2], list(reads))
        self.assertEqual([read1, read2], list(reads))

    def testSubclassWithAdditionalReads(self):
        """
        A Reads subclass with an iter method that is then added to manually
        must result in an instance with a correct iterator.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        read3 = Read('id3', 'AC')

        class ReadsSubclass(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = ReadsSubclass()
        reads.add(read3)
        self.assertEqual(sorted([read1, read2, read3]), sorted(reads))

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
        error = "Format must be either 'fasta', 'fastq' or 'fasta-ss'\\."
        six.assertRaisesRegex(self, ValueError, error, reads.save, 'file',
                              'xxx')
        # The output file must not exist following the save() failure.
        error = "No such file or directory: 'file'"
        six.assertRaisesRegex(self, OSError, error, stat, 'file')

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

    def testSaveReturnsReadCount(self):
        """
        The save method on a Reads instance must return the number
        of reads in the instance.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            result = reads.save('filename')
            self.assertIs(2, result)

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
        six.assertRaisesRegex(self, ValueError, error, reads.save, 'file',
                              'fastq')
        # The output file must not exist following the save() failure.
        error = "No such file or directory: 'file'"
        six.assertRaisesRegex(self, OSError, error, stat, 'file')

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

    def testUnfilteredLengthBeforeIterating(self):
        """
        A Reads instance must raise RuntimeError if its unfilteredLength method
        is called before it has been iterated.
        """
        reads = Reads()
        error = ('^The unfiltered length of a Reads instance is unknown until '
                 'it has been iterated\\.$')
        six.assertRaisesRegex(self, RuntimeError, error,
                              reads.unfilteredLength)

    def testUnfilteredLengthNoReads(self):
        """
        A Reads instance with no reads must have an unfiltered length of zero.
        """
        reads = Reads()
        list(reads)
        self.assertEqual(0, reads.unfilteredLength())

    def testUnfilteredLengthAdditionalReads(self):
        """
        A Reads instance that has been added to manually must have the correct
        unfiltered length.
        """
        reads = Reads()
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        list(reads)
        self.assertEqual(2, reads.unfilteredLength())

    def testUnfilteredLengthAdditionalReadsAfterFiltering(self):
        """
        A Reads instance that has been added to manually and then filtered must
        have the correct (original) unfiltered length and the filtered list it
        returns must be correct.
        """
        reads = Reads()
        read1 = Read('id1', 'ATTA')
        read2 = Read('id2', 'AC')
        reads.add(read1)
        reads.add(read2)
        reads.filter(minLength=3)
        self.assertEqual([read1], list(reads))
        self.assertEqual(2, reads.unfilteredLength())

    def testUnfilteredLengthInitialReads(self):
        """
        A Reads instance that has been given reads initially must have the
        correct unfiltered length.
        """
        read1 = Read('id1', 'AT')
        read2 = Read('id2', 'AC')
        reads = Reads([read1, read2])
        list(reads)
        self.assertEqual(2, reads.unfilteredLength())

    def testUnfilteredLengthInitialReadsAfterFiltering(self):
        """
        A Reads instance that has been given reads intially and then filtered
        must have the correct (original) unfiltered length and the filtered
        list it returns must be correct.
        """
        read1 = Read('id1', 'ATTA')
        read2 = Read('id2', 'AC')
        reads = Reads([read1, read2])
        reads.filter(minLength=3)
        self.assertEqual([read1], list(reads))
        self.assertEqual(2, reads.unfilteredLength())

    def testUnfilteredLengthInitialReadsIsReads(self):
        """
        A Reads instance that has been given another filtered Reads instance
        intially must have the correct (original) unfiltered length and the
        filtered list it returns must be correct.
        """
        read1 = Read('id1', 'ATTA')
        read2 = Read('id2', 'AC')
        initialReads = Reads([read1, read2])
        initialReads.filter(minLength=3)

        reads = Reads(initialReads)
        self.assertEqual([read1], list(reads))
        self.assertEqual(2, reads.unfilteredLength())

    def testUnfilteredLengthInitialReadsIsReadsWithAdditional(self):
        """
        A Reads instance that has been given another filtered Reads instance
        intially and then an additional read must have the correct (original)
        unfiltered length and the filtered list it returns must be correct.
        """
        read1 = Read('id1', 'ATTA')
        read2 = Read('id2', 'AC')
        initialReads = Reads([read1, read2])
        initialReads.filter(minLength=3)

        read3 = Read('id3', 'AC')
        reads = Reads(initialReads)
        reads.add(read3)
        self.assertEqual(sorted((read1, read3)), sorted(reads))
        self.assertEqual(3, reads.unfilteredLength())

    def testUnfilteredLengthInitialSubclassWithNoLen(self):
        """
        If a Reads instance is given a Reads subclass (with no __len__)
        instance intially, it must have the correct (original)
        unfiltered length and the filtered list it returns must be correct.
        """
        read1 = Read('id1', 'ATTA')
        read2 = Read('id2', 'AC')

        class Subclass(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = Reads(Subclass())
        self.assertEqual(sorted((read1, read2)), sorted(reads))
        self.assertEqual(2, reads.unfilteredLength())

    def testUnfilteredLengthInitialSubclassThenFiltered(self):
        """
        If a Reads instance is given a Reads subclass instance intially and is
        then filtered, it must have the correct (original) unfiltered length
        and the filtered list it returns must be correct.
        """
        read1 = Read('id1', 'ATTA')
        read2 = Read('id2', 'AC')

        class Subclass(Reads):
            def iter(self):
                yield read1
                yield read2

        reads = Reads(Subclass())
        reads.filter(maxLength=3)
        self.assertEqual([read2], sorted(reads))
        self.assertEqual(2, reads.unfilteredLength())

    def testUnfilteredLengthInitialSubclassWithAdditionalThenFiltered(self):
        """
        If a Reads instance is given a Reads subclass instance that has been
        added to intially and is then filtered, it must have the correct
        (original) unfiltered length and the filtered list it returns must be
        correct.
        """
        read1 = Read('id1', 'ATTA')
        read2 = Read('id2', 'AC')
        read3 = Read('id3', 'AC')

        class Subclass(Reads):
            def iter(self):
                yield read1
                yield read2

        initial = Subclass()
        initial.add(read3)
        reads = Reads(initial)
        reads.filter(maxLength=3)
        self.assertEqual(sorted([read2, read3]), sorted(reads))
        self.assertEqual(3, reads.unfilteredLength())


class TestReadsFiltering(TestCase):
    """
    Tests of filtering dark.reads.Reads instances.

    These tests use the convenience 'filter' method on the Reads class, which
    makes a ReadFilter instance. For that reason there are no explicit separate
    tests of the C{ReadFilter} class (there probably should be, but all the
    filtering tests below were written before the filtering code in C{Reads}
    was pulled out into a separate C{ReadFilter} class. All forms of filtering
    are tested.
    """

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

    def testFilterReturnsReadInstance(self):
        """
        The filter method must return a C{Reads} instance.
        """
        self.assertTrue(isinstance(Reads().filter(), Reads))

    def testFilteredReadsInstanceHasExpectedLength(self):
        """
        After filtering, the returned Reads instance must have the expected
        length.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ACG')
        read3 = Read('id3', 'AC')
        read4 = Read('id4', 'A')
        reads.add(read1)
        reads.add(read2)
        reads.add(read3)
        reads.add(read4)
        result = reads.filter(minLength=3)
        self.assertEqual(2, len(list(result)))

    def testAddFiltersThenClearFilters(self):
        """
        If filters are added and then all filters are cleared, the result must
        be the same as the reads that were originally added.
        """
        initial = [Read('id1', 'ATCG'), Read('id2', 'ACG'), Read('id3', 'AC'),
                   Read('id4', 'A')]
        reads = Reads(initial)
        result = reads.filter(minLength=3).filter(maxLength=3).clearFilters()
        self.assertEqual(initial, list(result))

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

    def testFilterRemoveGapsWithQuality(self):
        """
        Filtering must be able to remove gaps, treating the quality string
        properly.
        """
        reads = Reads()
        reads.add(Read('id', '-AT--CG-', '12345678'))
        result = reads.filter(removeGaps=True)
        self.assertEqual([Read('id', 'ATCG', '2367')], list(result))

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

    def testFilterKeepSequencesNoSequences(self):
        """
        Filtering must be able to filter reads based on their sequential
        number when no sequences are wanted.
        """
        reads = Reads((Read('cow', 'A'), Read('dog', 'G'), Read('cat', 'T')))
        result = reads.filter(keepSequences=set())
        self.assertEqual([], list(result))

    def testFilterKeepSequences(self):
        """
        Filtering must be able to filter reads based on their sequential
        number.
        """
        reads = Reads()
        reads.add(Read('cow', 'AA'))
        reads.add(Read('dog', 'GG'))
        reads.add(Read('cat', 'TT'))
        result = reads.filter(keepSequences=set([0, 2]))
        self.assertEqual([Read('cow', 'AA'), Read('cat', 'TT')], list(result))

    def testFilterRemoveSequencesNoSequences(self):
        """
        Filtering must be able to filter reads based on their sequential
        number when no sequences are excluded.
        """
        animals = [Read('cow', 'A'), Read('dog', 'G'), Read('cat', 'T')]
        reads = Reads(animals)
        result = reads.filter(removeSequences=set())
        self.assertEqual(animals, list(result))

    def testFilterRemoveSequences(self):
        """
        Filtering must be able to exclude reads based on their sequential
        number.
        """
        reads = Reads()
        reads.add(Read('cow', 'AA'))
        reads.add(Read('dog', 'GG'))
        reads.add(Read('cat', 'TT'))
        result = reads.filter(removeSequences=set([1]))
        self.assertEqual([Read('cow', 'AA'), Read('cat', 'TT')], list(result))

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
        Filtering on sequence duplicates must work correctly. The first of a
        set of duplicated reads is the one that should be retained.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(removeDuplicates=True)
        self.assertEqual([read1], list(result))

    def testFilterDuplicatesById(self):
        """
        Filtering on read id duplicates must work correctly. The first of a
        set of duplicated reads is the one that should be retained.
        """
        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id1', 'ATTT')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(removeDuplicatesById=True)
        self.assertEqual([read1], list(result))

    def testFilterRemoveDescriptions(self):
        """
        Removing read id descriptions must work correctly.
        """
        read1 = Read('id1 description', 'ATCG')
        read2 = Read('id2 another description', 'ATTT')
        read3 = Read('id3', 'ATT')
        reads = Reads([read1, read2, read3])
        result = reads.filter(removeDescriptions=True)
        self.assertEqual(['id1', 'id2', 'id3'],
                         [read.id for read in result])

    def testFilterDoNotRemoveDescriptions(self):
        """
        Not removing read id descriptions must work correctly.
        """
        read1 = Read('id1 description', 'ATCG')
        read2 = Read('id2 another description', 'ATTT')
        read3 = Read('id3', 'ATT')
        reads = Reads([read1, read2, read3])
        result = reads.filter(removeDescriptions=False)
        self.assertEqual([read1, read2, read3], list(result))

    def testFilterWithModifierThatOmits(self):
        """
        Filtering with a modifier function must work correctly if the modifier
        returns C{None} for some reads.
        """
        def modifier(read):
            if read.id == 'id1':
                return read

        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(modifier=modifier)
        self.assertEqual([read1], list(result))

    def testFilterWithModifierThatChangesIds(self):
        """
        Filtering with a modifier function must work correctly if the modifier
        changes the ids of reads.
        """
        def modifier(read):
            read.id = read.id.upper()
            return read

        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(modifier=modifier)
        self.assertEqual([Read('ID1', 'ATCG'), Read('ID2', 'ATCG')],
                         list(result))

    def testFilterWithModifierThatOmitsAndChangesIds(self):
        """
        Filtering with a modifier function must work correctly if the modifier
        omits some reads and changes the ids of others.
        """
        def modifier(read):
            if read.id == 'id1':
                read.id = 'xxx'
                return read

        reads = Reads()
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads.add(read1)
        reads.add(read2)
        result = reads.filter(modifier=modifier)
        self.assertEqual([Read('xxx', 'ATCG')], list(result))

    def testFilterRandomSubsetSizeZeroNoReads(self):
        """
        Asking for a random subset of length zero must work as expected when
        there are no reads in the Reads instance.
        """
        self.assertEqual([],
                         list(Reads().filter(randomSubset=0, trueLength=0)))

    def testFilterRandomSubsetSizeZeroTwoReads(self):
        """
        Asking for a random subset of length zero from a set of two reads must
        work as expected.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads = Reads(initialReads=[read1, read2])
        result = reads.filter(randomSubset=0, trueLength=2)
        self.assertEqual([], list(result))

    def testFilterRandomSubsetOfZeroReads(self):
        """
        Asking for a non-zero random subset of a set of zero reads must work
        as expected.
        """
        reads = Reads()
        result = reads.filter(randomSubset=5, trueLength=0)
        self.assertEqual([], list(result))

    def testFilterRandomSubsetOfOneFromOneRead(self):
        """
        Asking for a size one random subset of a set of one read must work
        as expected.
        """
        read = Read('id', 'ATCG')
        reads = Reads(initialReads=[read])
        result = reads.filter(randomSubset=1, trueLength=1)
        self.assertEqual([read], list(result))

    def testFilterRandomSubsetOfFiveFromOneRead(self):
        """
        Asking for a size five random subset of a set of one read must work
        as expected.
        """
        read = Read('id', 'ATCG')
        reads = Reads(initialReads=[read])
        result = reads.filter(randomSubset=5, trueLength=1)
        self.assertEqual([read], list(result))

    def testFilterRandomSubsetOfFiveFromFiveReads(self):
        """
        Asking for a size five random subset of a set of five reads must work
        as expected.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATCG')
        read4 = Read('id4', 'ATCG')
        read5 = Read('id5', 'ATCG')
        reads = Reads(initialReads=[read1, read2, read3, read4, read5])
        result = reads.filter(randomSubset=5, trueLength=5)
        self.assertEqual([read1, read2, read3, read4, read5], list(result))

    def testFilterRandomSubsetOfTwoFromFiveReads(self):
        """
        Asking for a size two random subset of a set of five reads must return
        two (different) reads.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATCG')
        read4 = Read('id4', 'ATCG')
        read5 = Read('id5', 'ATCG')
        reads = Reads(initialReads=[read1, read2, read3, read4, read5])
        result = reads.filter(randomSubset=2, trueLength=5)
        self.assertEqual(2, len(set(result)))

    def testSampleFractionAndRandomSubsetRaisesValueError(self):
        """
        Asking for filtering of a sample fraction and a random subset at the
        same time must raise a ValueError.
        """
        reads = Reads()
        error = ("^randomSubset and sampleFraction cannot be used "
                 "simultaneously in a filter. Make two read filters "
                 "instead\\.$")
        six.assertRaisesRegex(self, ValueError, error, reads.filter,
                              sampleFraction=0.1, randomSubset=3)

    def testSampleFractionAndNoTrueLengthRaisesValueError(self):
        """
        Asking for filtering of a sample fraction without passing a trueLength
        must raise a ValueError.
        """
        reads = Reads()
        error = "^trueLength must be supplied if randomSubset is specified\\.$"
        six.assertRaisesRegex(self, ValueError, error, reads.filter,
                              randomSubset=3)

    def testSampleFractionZero(self):
        """
        Asking for a sample fraction of 0.0 from a set of five reads must
        return the empty list.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATCG')
        read4 = Read('id4', 'ATCG')
        read5 = Read('id5', 'ATCG')
        reads = Reads(initialReads=[read1, read2, read3, read4, read5])
        result = reads.filter(sampleFraction=0.0)
        self.assertEqual(0, len(list(result)))

    def testSampleFractionOne(self):
        """
        Asking for a sample fraction of 1.0 from a set of five reads must
        return all reads.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATCG')
        read4 = Read('id4', 'ATCG')
        read5 = Read('id5', 'ATCG')
        reads = Reads(initialReads=[read1, read2, read3, read4, read5])
        result = reads.filter(sampleFraction=1.0)
        self.assertEqual([read1, read2, read3, read4, read5], list(result))

    def testSampleFractionPointOne(self):
        """
        Asking for a sample fraction of 0.1 from a set of 100 reads must
        return 11 reads (given a particular random seed value).
        """
        seed(1)
        reads = Reads(initialReads=[Read('id1', 'ATCG')] * 100)
        result = reads.filter(sampleFraction=0.1)
        self.assertEqual(11, len(list(result)))

    def testLineNumberFileFirstLineTooSmall(self):
        """
        If a line number file is passed to filter but its first line is
        less than 1, a ValueError must be raised.
        """
        data = '0\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads([])
            error = ("^First line of sequence number file 'file' must be at "
                     "least 1\\.$")
            with six.assertRaisesRegex(self, ValueError, error):
                reads.filter(sequenceNumbersFile='file')

    def testLineNumberFileNonAscending(self):
        """
        If a line number file is passed to filter but it contains non-ascending
        line numbers, a ValueError must be raised.
        """
        data = '2\n2\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            read1 = Read('id1', 'ATCG')
            read2 = Read('id2', 'ATCG')
            read3 = Read('id3', 'ATCG')
            reads = Reads(initialReads=[read1, read2, read3])
            error = ("^Line number file 'file' contains non-ascending numbers "
                     "2 and 2\\.$")
            with six.assertRaisesRegex(self, ValueError, error):
                list(reads.filter(sequenceNumbersFile='file'))

    def testLineNumberFileEmpty(self):
        """
        If an empty line number file is passed to filter, no sequences should
        be returned.
        """
        data = ''
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            read1 = Read('id1', 'ATCG')
            read2 = Read('id2', 'ATCG')
            reads = Reads(initialReads=[read1, read2])
            result = reads.filter(sequenceNumbersFile='file')
            self.assertEqual([], list(result))

    def testLineNumberFile(self):
        """
        If a line number file is passed to filter, the correct sequences should
        be returned.
        """
        data = '1\n3\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            read1 = Read('id1', 'ATCG')
            read2 = Read('id2', 'ATCG')
            read3 = Read('id3', 'ATCG')
            read4 = Read('id4', 'ATCG')
            reads = Reads(initialReads=[read1, read2, read3, read4])
            result = reads.filter(sequenceNumbersFile='file')
            self.assertEqual([read1, read3], list(result))

    def testLineNumberFileRunOutOfSequences(self):
        """
        If a line number file is passed to filter, and it contains numbers that
        are bigger than the number of sequences, the correct sequences should
        be returned.
        """
        data = '1\n3\n200\n'
        mockOpener = mockOpen(read_data=data)
        with patch.object(builtins, 'open', mockOpener):
            read1 = Read('id1', 'ATCG')
            read2 = Read('id2', 'ATCG')
            read3 = Read('id3', 'ATCG')
            read4 = Read('id4', 'ATCG')
            reads = Reads(initialReads=[read1, read2, read3, read4])
            result = reads.filter(sequenceNumbersFile='file')
            self.assertEqual([read1, read3], list(result))

    def testKeepSites(self):
        """
        If only a certain set of sites should be kept, the correct sequences
        should be returned.
        """
        read1 = Read('id1', 'ATCGAT')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATC')
        read4 = Read('id4', 'AT')
        reads = Reads(initialReads=[read1, read2, read3, read4])
        result = reads.filter(keepSites={1, 2})
        self.assertEqual(
            [
                Read('id1', 'TC'),
                Read('id2', 'TC'),
                Read('id3', 'TC'),
                Read('id4', 'T'),
            ],
            list(result))

    def testKeepSitesNoSites(self):
        """
        If the empty set of sites should be kept, the correct (empty)
        sequences should be returned.
        """
        read1 = Read('id1', 'ATCGAT')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATC')
        read4 = Read('id4', 'AT')
        reads = Reads(initialReads=[read1, read2, read3, read4])
        result = reads.filter(keepSites=set())
        self.assertEqual(
            [
                Read('id1', ''),
                Read('id2', ''),
                Read('id3', ''),
                Read('id4', ''),
            ],
            list(result))

    def testKeepSitesAllSites(self):
        """
        If the set of all sites should be kept, the correct (full) sequences
        should be returned.
        """
        read1 = Read('id1', 'ATCGAT')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATC')
        read4 = Read('id4', 'AT')
        reads = Reads(initialReads=[read1, read2, read3, read4])
        result = reads.filter(keepSites=set(range(6)))
        self.assertEqual([read1, read2, read3, read4], list(result))

    def testKeepSitesWithQuality(self):
        """
        If only a certain set of sites should be kept, the correct sequences
        should be returned, including a modified quality string.
        """
        reads = Reads(initialReads=[Read('id1', 'ATCGAT', '123456')])
        result = reads.filter(keepSites={1, 2})
        self.assertEqual([Read('id1', 'TC', '23')], list(result))

    def testKeepSitesOutOfRange(self):
        """
        If only a certain set of sites should be kept, but the kept sites
        are higher than the length of the input sequences, the correct
        (empty) sequences should be returned.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCGCC')
        reads = Reads(initialReads=[read1, read2])
        result = reads.filter(keepSites={100, 200})
        self.assertEqual(
            [
                Read('id1', ''),
                Read('id2', ''),
            ],
            list(result))

    def testRemoveSites(self):
        """
        If a certain set of sites should be removed, the correct sequences
        should be returned.
        """
        read1 = Read('id1', 'ATCGCC')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATC')
        read4 = Read('id4', 'A')
        reads = Reads(initialReads=[read1, read2, read3, read4])
        result = reads.filter(removeSites={1, 2})
        self.assertEqual(
            [
                Read('id1', 'AGCC'),
                Read('id2', 'AG'),
                Read('id3', 'A'),
                Read('id4', 'A'),
            ],
            list(result))

    def testRemoveSitesNoSites(self):
        """
        If the empty set of sites should be removed, the correct (full)
        sequences should be returned.
        """
        read1 = Read('id1', 'ATCGCC')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATC')
        read4 = Read('id4', 'A')
        reads = Reads(initialReads=[read1, read2, read3, read4])
        result = reads.filter(removeSites=set())
        self.assertEqual([read1, read2, read3, read4],
                         list(result))

    def testRemoveSitesAllSites(self):
        """
        If the full set of sites should be removed, the correct (empty)
        sequences should be returned.
        """
        read1 = Read('id1', 'ATCGCC')
        read2 = Read('id2', 'ATCG')
        read3 = Read('id3', 'ATC')
        read4 = Read('id4', 'A')
        reads = Reads(initialReads=[read1, read2, read3, read4])
        result = reads.filter(removeSites=set(range(6)))
        self.assertEqual(
            [
                Read('id1', ''),
                Read('id2', ''),
                Read('id3', ''),
                Read('id4', ''),
            ],
            list(result))

    def testRemoveSitesWithQuality(self):
        """
        If a certain set of sites should be removed, the correct sequences
        should be returned, including a modified quality string.
        """
        reads = Reads(initialReads=[Read('id1', 'ATCGAT', '123456')])
        result = reads.filter(removeSites={1, 2})
        self.assertEqual([Read('id1', 'AGAT', '1456')], list(result))

    def testRemoveSitesOutOfRange(self):
        """
        If a certain set of sites should be removed, but the sites
        are higher than the length of the input sequences, the correct
        (full) sequences should be returned.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCGAA')
        reads = Reads(initialReads=[read1, read2])
        result = reads.filter(removeSites={100, 200})
        self.assertEqual([read1, read2], list(result))

    def testRemoveAndKeepSites(self):
        """
        Passing a keepSites and a removeSites set to reads.filter
        must result in a ValueError.
        """
        error = ('^Cannot simultaneously filter using keepSites and '
                 'removeSites\\. Call filter twice in succession instead\\.$')
        six.assertRaisesRegex(self, ValueError, error, Reads().filter,
                              keepSites={4}, removeSites={5})

    def testIdLambda(self):
        """
        A passed idLambda function should produce the expected read ids.
        """
        read = Read('id1', 'ATCGCC')
        reads = Reads(initialReads=[read])
        result = reads.filter(idLambda='lambda id: "x-" + id.upper()')
        self.assertEqual('x-ID1', list(result)[0].id)

    def testIdLambdaReturningNone(self):
        """
        A passed idLambda function should produce the expected read ids,
        including when it returns None.
        """
        read1 = Read('id1', 'ATCGCC')
        read2 = Read('id2', 'GGATCG')
        reads = Reads(initialReads=[read1, read2])
        result = reads.filter(
            idLambda='lambda id: "aa" if id.find("1") > -1 else None')
        (result,) = list(result)
        self.assertEqual('aa', result.id)

    def testReadLambda(self):
        """
        A passed readLambda function should produce the expected reads.
        """
        read = Read('id1', 'ATCGCC')
        reads = Reads(initialReads=[read])
        result = reads.filter(readLambda='lambda r: Read("hey", "AAA")')
        (result,) = list(result)
        self.assertEqual(Read('hey', 'AAA'), result)

    def testReadLambdaReturningNone(self):
        """
        A passed readLambda function should produce the expected reads,
        including when it returns None.
        """
        read1 = Read('xid1', 'ATCGCC')
        read2 = Read('yid2', 'GGATCG')
        reads = Reads(initialReads=[read1, read2])
        result = reads.filter(
            readLambda=('lambda r: Read(r.id + "-x", r.sequence[:2]) '
                        'if r.id.startswith("x") else None'))
        (result,) = list(result)
        self.assertEqual(Read('xid1-x', 'AT'), result)

    def testReverse(self):
        """
        When reverse=True, reads with the expected sequences and qualities
        must be returned.
        """
        read1 = Read('id1', 'ATCGCC', '123456')
        read2 = Read('id2', 'GGATCG', '987654')
        reads = Reads(initialReads=[read1, read2])
        result1, result2 = list(reads.filter(reverse=True))
        self.assertEqual(Read('id1', 'CCGCTA', '654321'), result1)
        self.assertEqual(Read('id2', 'GCTAGG', '456789'), result2)

    def testReverseComplement(self):
        """
        When reverseComplement=True, reads with the expected sequences and
        qualities must be returned when the reads are of type DNARead.
        """
        read1 = DNARead('id1', 'ATCGCC', '123456')
        read2 = DNARead('id2', 'GGATCG', '987654')
        reads = Reads(initialReads=[read1, read2])
        result1, result2 = list(reads.filter(reverseComplement=True))
        self.assertEqual(Read('id1', 'GGCGAT', '654321'), result1)
        self.assertEqual(Read('id2', 'CGATCC', '456789'), result2)

    def testReverseAndReverseComplement(self):
        """
        When reverseComplement=True and reverse=True, the expected reads must
        be returned (reverse complemented) when the reads are of type DNARead.
        """
        read1 = DNARead('id1', 'ATCGCC', '123456')
        read2 = DNARead('id2', 'GGATCG', '987654')
        reads = Reads(initialReads=[read1, read2])
        result1, result2 = list(reads.filter(reverseComplement=True,
                                             reverse=True))
        self.assertEqual(Read('id1', 'GGCGAT', '654321'), result1)
        self.assertEqual(Read('id2', 'CGATCC', '456789'), result2)

    def testReverseComplementNonDNA(self):
        """
        When reverseComplement=True and the read is not a DNARead, an
        AttributeError must be raised.
        """
        reads = Reads(initialReads=[Read('id1', 'ATCGCC', '123456')])
        error = "^'Read' object has no attribute 'reverseComplement'$"
        six.assertRaisesRegex(self, AttributeError, error, list,
                              reads.filter(reverseComplement=True))

    def testReverseComplementAARead(self):
        """
        When reverseComplement=True and the read is an AARead, an
        AttributeError must be raised.
        """
        reads = Reads(initialReads=[AARead('id1', 'ATCGCC', '123456')])
        error = "^'AARead' object has no attribute 'reverseComplement'$"
        six.assertRaisesRegex(self, AttributeError, error, list,
                              reads.filter(reverseComplement=True))


class TestReadsInRAM(TestCase):
    """
    Test the ReadsInRAM class.
    """

    def testNoReads(self):
        """
        A ReadsInRAM instance with no reads must return an empty iterator.
        """
        reads = ReadsInRAM()
        self.assertEqual([], list(reads))

    def testAdd(self):
        """
        It must be possible to add reads to a ReadsInRAM instance.
        """
        reads = ReadsInRAM()
        read = Read('id', 'ACGT')
        reads.add(read)
        self.assertEqual([read], list(reads))

    def testOneReadLength(self):
        """
        A ReadsInRAM instance with one read must have length one.
        """
        read1 = Read('id1', 'ATCG')
        reads = ReadsInRAM([read1])
        self.assertEqual(1, len(reads))

    def testOneReadList(self):
        """
        A ReadsInRAM instance with one read must iterate as expected.
        """
        read1 = Read('id1', 'ATCG')
        reads = ReadsInRAM([read1])
        self.assertEqual([read1], list(reads))

    def testOneReadIndex(self):
        """
        A ReadsInRAM instance with one read must be able to be indexed.
        """
        read1 = Read('id1', 'ATCG')
        reads = ReadsInRAM([read1])
        self.assertEqual(read1, reads[0])

    def testTwoReadsLength(self):
        """
        A ReadsInRAM instance with two reads must have length one.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads = ReadsInRAM([read1, read2])
        self.assertEqual(2, len(reads))

    def testTwoReadsList(self):
        """
        A ReadsInRAM instance with two reads must iterate as expected.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads = ReadsInRAM([read1, read2])
        self.assertEqual([read1, read2], list(reads))

    def testTwoReadsIndex(self):
        """
        A ReadsInRAM instance with two reads must be able to be indexed.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads = ReadsInRAM([read1, read2])
        self.assertEqual(read1, reads[0])
        self.assertEqual(read2, reads[1])

    def testFromReads(self):
        """
        A ReadsInRAM instance must be able to initialize itself from a Reads
        instance.
        """
        read1 = Read('id1', 'ATCG')
        read2 = Read('id2', 'ATCG')
        reads = Reads([read1, read2])
        readsInRAM = ReadsInRAM(reads)
        self.assertEqual(2, len(readsInRAM))
        self.assertEqual(read1, readsInRAM[0])
        self.assertEqual(read2, readsInRAM[1])

    def testFastaFile(self):
        """
        A ReadsInRAM instance should be able to be initialized from a
        FastaReads instance and be iterated twice without the underlying file
        being opened twice.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('file1.fasta', filename)
                    self.count += 1
                    return StringIO('>id1\nACTG\n>id2\nAA\n')
                else:
                    self.test.fail('We are only supposed to be called once!')

        sideEffect = SideEffect(self)
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = FastaReads('file1.fasta')
            readsInRAM = ReadsInRAM(reads)
            expected = [Read('id1', 'ACTG'), Read('id2', 'AA')]
            self.assertEqual(expected, list(readsInRAM))
            self.assertEqual(expected, list(readsInRAM))
            self.assertEqual(1, sideEffect.count)

    def testSetItem(self):
        """
        It must be possible to set a value for a ReadsInRAM index.
        """
        read1 = Read('id1', 'ATCG')
        reads = ReadsInRAM([read1])
        read2 = Read('id2', 'ATCG')
        reads[0] = read2
        self.assertEqual(read2, reads[0])


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


class TestSitesMatching(TestCase):
    """
    Tests for the Reads.sitesMatching method.
    """
    def testNoMatches(self):
        """
        If no reads match the target bases, sitesMatching must return the
        empty set.
        """
        reads = Reads()
        reads.add(Read('id1', 'aaaaaa'))
        result = reads.sitesMatching({'b'}, matchCase=True, any_=True)
        self.assertEqual(set(), result)

    def testAllMatches(self):
        """
        If a read's sites all match the target bases, sitesMatching must
        return the full set of sites.
        """
        reads = Reads()
        reads.add(Read('id1', 'aaaga'))
        result = reads.sitesMatching({'a', 'g'}, matchCase=True, any_=True)
        self.assertEqual(set(range(5)), result)

    def testPartialMatch(self):
        """
        If some of a read's sites all match the target bases, sitesMatching
        must return the expected set of sites.
        """
        reads = Reads()
        reads.add(Read('id1', 'aaaga'))
        result = reads.sitesMatching({'a'}, matchCase=True, any_=True)
        self.assertEqual({0, 1, 2, 4}, result)

    def testMatchCase(self):
        """
        If only some of a read's sites all match the target bases with
        matching case, sitesMatching must return the expected set of sites.
        """
        reads = Reads()
        reads.add(Read('id1', 'aAAga'))
        result = reads.sitesMatching({'a'}, matchCase=True, any_=False)
        self.assertEqual({0, 4}, result)

    def testIgnoreCase(self):
        """
        If only some of a read's sites all match the target bases with
        matching case, sitesMatching must return the expected set of sites
        when we tell it to ignore case.
        """
        reads = Reads()
        reads.add(Read('id1', 'aAAga'))
        result = reads.sitesMatching({'a'}, matchCase=False, any_=False)
        self.assertEqual({0, 1, 2, 4}, result)

    def testMultipleReadsAny(self):
        """
        If multiple reads are given, sitesMatching must return the expected
        set of sites when we pass any_=True.
        """
        reads = Reads()
        reads.add(Read('id1', 'aaaga'))
        reads.add(Read('id2', 'caata'))
        result = reads.sitesMatching({'c', 'g'}, matchCase=False, any_=True)
        self.assertEqual({0, 3}, result)

    def testMultipleReadsAll(self):
        """
        If multiple reads are given, sitesMatching must return the expected
        set of sites when we pass any_=False.
        """
        reads = Reads()
        reads.add(Read('id1', 'aaaga'))
        reads.add(Read('id2', 'caaga'))
        result = reads.sitesMatching({'c', 'g'}, matchCase=False, any_=False)
        self.assertEqual({3}, result)

    def testMultipleReadsAllWithDifferingLengths(self):
        """
        If multiple reads are given, sitesMatching must return the expected
        set of sites when we pass any_=False and the reads have differing
        lengths.
        """
        reads = Reads()
        reads.add(Read('id2', 'caa-agt-'))
        reads.add(Read('id1', '-aa-a'))
        result = reads.sitesMatching({'-'}, matchCase=False, any_=False)
        self.assertEqual({3}, result)

    def testMultipleReadsAnyWithDifferingLengths(self):
        """
        If multiple reads are given, sitesMatching must return the expected
        set of sites when we pass any_=True and the reads have differing
        lengths.
        """
        reads = Reads()
        reads.add(Read('id2', 'caa-agt-'))
        reads.add(Read('id1', '-aa-a'))
        result = reads.sitesMatching({'-'}, matchCase=False, any_=True)
        self.assertEqual({0, 3, 7}, result)

from unittest import TestCase, skipUnless
from contextlib import contextmanager
from pathlib import Path
from tempfile import mkdtemp

from dark.cigar import makeCigar
from dark.consensus import consensusFromBAM
from dark.process import Executor
from dark.reads import DNARead
from dark.sam import (
    samtoolsInstalled, UnequalReferenceLengthError, UnknownReference,
    UnspecifiedReference)
from dark.utils import matchOffset

# From https://samtools.github.io/hts-specs/SAMv1.pdf
CINS, CDEL, CMATCH = 'IDM'


@contextmanager
def makeBAM(template, secondReference=None):
    """
    A context manager decorator to make a simple BAM file from a template.
    Note that this code invokes samtools.

    @param template: An iterable of C{str} sequences. The first will be treated
        as the reference, and then subsequent pairs (if any) will be treated as
        read and quality strings. Reads and quality strings can be indented
        with spaces to show where the read aligns with the reference.
    @param secondReference: If not C{None}, the C{str} id of a second reference
        to add to the BAM header.
    @return: A context manager that produces a 2-tuple containing the reference
        C{DNARead} instance and the C{Path} of the BAM file.
    """
    if len(template) % 2 != 1:
        raise ValueError(
            'The template must have an odd number of strings, specifying the '
            'reference sequence, then zero or more read/quality pairs.')

    refId = 'ref-id'
    leftPaddedReference = template[0]
    reference = DNARead(refId, leftPaddedReference.lstrip().replace('-', ''))
    nSeqs = (len(template) - 1) >> 1
    dirname = mkdtemp(prefix='test-consensus-')
    e = Executor()

    try:
        samFile = Path(dirname) / 'file.sam'
        bamFile = Path(dirname) / 'file.bam'
        with open(samFile, 'w') as fp:
            print(f'@SQ\tSN:{refId}\tLN:{len(reference)}', file=fp)
            if secondReference:
                # Add a second reference, of an arbitrary length.
                print(f'@SQ\tSN:{secondReference}\tLN:100', file=fp)

            for count in range(nSeqs):
                leftPaddedQuery = template[count * 2 + 1].rstrip()
                leftPaddedQuality = template[count * 2 + 2].rstrip()
                assert len(leftPaddedQuery) == len(leftPaddedQuality)
                query = leftPaddedQuery.lstrip()
                quality = leftPaddedQuality.lstrip()
                queryNoGaps = qualityNoGaps = ''
                for queryBase, qualityBase in zip(query, quality):
                    if queryBase != '-':
                        queryNoGaps += queryBase
                        qualityNoGaps += qualityBase

                print('\t'.join(map(str, (
                    f'read{count}',  # QNAME (query name)
                    0,  # FLAGS
                    refId,  # RNAME (reference name)
                    matchOffset(leftPaddedReference, leftPaddedQuery) + 1,
                    30,  # MAPQ (mapping quality)
                    makeCigar(leftPaddedReference, leftPaddedQuery),  # CIGAR
                    '*',  # MRNM (mate reference name)
                    0,  # MPOS (mate position)
                    0,  # ISIZE (insert size)
                    queryNoGaps,  # SEQ
                    qualityNoGaps,  # QUAL
                ))), file=fp)

        e.execute(f'samtools view -b -o {str(bamFile)!r} {str(samFile)!r}')
        e.execute(f'samtools index {str(bamFile)!r}')
        yield (reference, bamFile)
    finally:
        e.execute(f'rm -fr {dirname!r}')


@skipUnless(samtoolsInstalled(), 'samtools is not installed')
class TestReferenceErrors(TestCase):
    """
    Test expected errors.
    """
    def testUnknownStrategy(self):
        """
        If an unknown strategy is passed, ValueError must be raised.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            strategy = 'xxx'
            error = rf'^Unknown consensus strategy {strategy!r}\.$'
            self.assertRaisesRegex(ValueError, error, consensusFromBAM,
                                   bamFilename, strategy=strategy)

    def testUnknownReferenceId(self):
        """
        If a reference id that is not in the BAM file is passed,
        UnknownReference must be raised.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            badId = 'bad-id'
            error = (
                rf'^BAM file {str(bamFilename)!r} does not mention a '
                rf'reference with id {badId!r}\. Known references are: '
                rf'ref-id\.$')
            self.assertRaisesRegex(UnknownReference, error, consensusFromBAM,
                                   bamFilename, referenceId=badId)

    def testUnknownReference(self):
        """
        If a reference sequence is passed and its id is not in the BAM
        file, UnknownReference must be raised.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            wrongReference = DNARead('unknown-id', template[0])
            error = (
                rf'^BAM file {str(bamFilename)!r} does not mention a '
                rf'reference with id {wrongReference.id!r}\. Known '
                rf'references are: ref-id\.$')
            self.assertRaisesRegex(UnknownReference, error, consensusFromBAM,
                                   bamFilename, reference=wrongReference)

    def testUnspecifiedReference(self):
        """
        If a reference is not specified and the BAM file mentions more than
        one, UnspecifiedReference must be raised.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template, secondReference='zzz') as (reference,
                                                          bamFilename):
            error = (
                rf'^BAM file {str(bamFilename)!r} mentions 2 references '
                rf'\(ref-id, zzz\) but you have not passed a referenceId '
                rf'argument or a reference sequence to indicate which one to '
                rf'use\.$')
            self.assertRaisesRegex(UnspecifiedReference, error,
                                   consensusFromBAM, bamFilename)

    def testReferenceOfWrongLength(self):
        """
        If a reference sequence is passed and it does not have the right
        length, UnequalReferenceLengthError must be raised.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            wrongReference = DNARead(reference.id, 'AA')
            error = (
                rf'^Reference with id {reference.id!r} has length 2, which '
                rf'does not match the length of reference {reference.id!r} '
                rf'\({len(template[0])}\) in BAM file '
                rf'{str(bamFilename)!r}\.$')
            self.assertRaisesRegex(UnequalReferenceLengthError, error,
                                   consensusFromBAM, bamFilename,
                                   referenceId=reference.id,
                                   reference=wrongReference)


class _Mixin:
    """
    Common (i.e., with and without considering FASTQ quality values) tests for
    consensuses making.
    """

    def testNoReadsReference(self):
        """
        If no reads are present and resolution of no-coverage bases is the
        reference sequence, the reference should be returned.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                template[0],
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 noCoverage='reference',
                                 ignoreQuality=self.ignoreQuality))

    def testNoReadsReferenceFromId(self):
        """
        The reference can be passed using just its id (as opposed to a full
        reference Read instance). A string of 'N's is returned because there
        are no reads and we pass noCoverage='N'.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'N' * len(template[0]),
                consensusFromBAM(bamFilename, referenceId=reference.id,
                                 noCoverage='N', lowCoverage='N',
                                 ignoreQuality=self.ignoreQuality))

    def testReferenceAndReferenceIdNotGiven(self):
        """
        If the reference to use is not given, the single reference in the BAM
        file should be used. A string of 'N's is returned because there
        are no reads and we pass noCoverage='N'.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'N' * len(template[0]),
                consensusFromBAM(bamFilename,
                                 noCoverage='N', lowCoverage='N',
                                 ignoreQuality=self.ignoreQuality))

    def testNoReadsN(self):
        """
        If no reads are present and resolution of no-coverage bases is 'N'
        (or '?', etc) a sequence of Ns (or ?s, etc) should be returned.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as (reference, bamFilename):
            for char in 'N?':
                self.assertEqual(
                    char * len(template[0]),
                    consensusFromBAM(bamFilename,
                                     reference=reference,
                                     noCoverage=char,
                                     ignoreQuality=self.ignoreQuality))

    def testLowReadsReference(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        bases is the reference sequence, the reference should be returned.
        """
        template = (
            'ACGTTCCG',
            '  GTT',
            '  ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                template[0],
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 lowCoverage='reference',
                                 minCoverage=2,
                                 ignoreQuality=self.ignoreQuality))

    def testLowReadsCharNoCoverageConsensus(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        bases is 'N' (or '? etc), then Ns (or ?s, etc) should be returned in
        the low coverage sites, and the reference in the sites with no
        coverage.
        """
        template = (
            'ACGTTCCG',
            '  GTT',
            '  ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            for char in 'N?':
                self.assertEqual(
                    'AC' + char * 3 + 'CCG',
                    consensusFromBAM(bamFilename,
                                     reference=reference,
                                     lowCoverage=char,
                                     minCoverage=2,
                                     ignoreQuality=self.ignoreQuality))

    def testLowReadsCharNoCoverageX(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        bases is 'N' (or '? etc), then Ns (or ?s, etc) should be returned in
        the low coverage sites, and (for example) 'X' in the sites with no
        coverage.
        """
        template = (
            'ACGTTCCG',
            '  GTT',
            '  ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            for char in 'N?':
                self.assertEqual(
                    'XX' + char * 3 + 'XXX',
                    consensusFromBAM(bamFilename,
                                     reference=reference,
                                     noCoverage='X',
                                     lowCoverage=char,
                                     minCoverage=2,
                                     ignoreQuality=self.ignoreQuality))

    def testOneReadMatchingPartOfTheReference(self):
        """
        If one read is present and it matches part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        reference should be returned.
        """
        template = (
            'ACGTTCCG',
            '  GTT',
            '  ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                template[0],
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 noCoverage='reference',
                                 ignoreQuality=self.ignoreQuality))

    def testOneReadDifferingFromPartOfTheReference(self):
        """
        If one read is present and it differs from part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        expected hybrid of the read and the reference should be returned.
        """
        template = (
            'ACGTTCCG',
            '  AAA',
            '  ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'ACAAACCG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 noCoverage='reference',
                                 ignoreQuality=self.ignoreQuality))

    def testTwoReadsDifferingFromPartOfTheReferenceSomeLowCoverage(self):
        """
        If two reads are present and they differ from part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        expected hybrid of the read and the reference should be returned, with
        the part of the reads that has insufficient coverage returning the
        low-coverage symbol (here '+').
        """
        template = (
            'ACGTTCCG',
            '  AAA',
            '  ???',
            '  AAAT',
            '  ????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'ACAAA+CG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 minCoverage=2,
                                 noCoverage='reference',
                                 lowCoverage='+',
                                 ignoreQuality=self.ignoreQuality))

    def testTwoReadsDifferingFromPartOfTheReferenceLowAndNoCoverage(self):
        """
        If two reads are present and they differ from part of the reference and
        resolution of no-coverage bases is 'N', the expected hybrid of the read
        and the reference should be returned, with the part of the reads that
        has insufficient coverage returning the low-coverage symbol (here '+').
        """
        template = (
            'ACGTTCCG',
            '  AAA',
            '  ???',
            '  AAAT',
            '  ????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'NNAAA+NN',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 minCoverage=2,
                                 noCoverage='N',
                                 lowCoverage='+',
                                 ignoreQuality=self.ignoreQuality))

    def testSimpleMajority(self):
        """
        If three reads result in a majority base at a site, that base should
        be in the consensus.
        """
        template = (
            'ACGT',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  C',
            '  ?',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'ACAT',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.5,
                                 ignoreQuality=self.ignoreQuality))

    def testSimpleMajorityBelowThreshold(self):
        """
        If conflicting reads at a site do not give a simple (above threshold)
        majority, the ambiguous code should be in the consensus.
        """
        template = (
            'ACGT',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  C',
            '  ?',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'ACMT',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testReadHasEarlierSites(self):
        """
        If a read has sites that come before the reference, they must
        be included in the consensus.
        """
        template = (
            ' ACGT',
            'AA',
            '??',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'AACGT',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testReadHasLaterSites(self):
        """
        If a read has sites that come after the reference, they must
        be included in the consensus.
        """
        template = (
            'ACGT',
            '   TAA',
            '   ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'ACGTAA',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testReadHasEarlierAndLaterSites(self):
        """
        If a read has sites that come before and after the reference, the sites
        must be included in the consensus.
        """
        template = (
            '  ACGT',
            'TCACGTGA',
            '????????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'TCACGTGA',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testReadsHaveEarlierAndLaterSites(self):
        """
        If two reads have sites that come before and after the reference, the
        sites must be included in the consensus. This is the same as the
        testReadHasEarlierAndLaterSites test above, but using two reads
        (instead of one) to give the before and after sites.
        """
        template = (
            '  ACGT',
            'TCA',
            '???',
            '     TGA',
            '     ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'TCACGTGA',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testOneDeletionFromReference(self):
        """
        A deletion from the reference must be handled correctly.
        """
        template = (
            'CACGTG',
            ' A-A',
            ' ?-?',
            ' A-A',
            ' ?-?',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'CAATG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testTwoDeletionsFromReference(self):
        """
        Two deletions from the reference must be handled correctly.
        """
        template = (
            'CACGTG',
            ' A-A-G',
            ' ?-?-?',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'CAAG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testOneInsertionInReference(self):
        """
        An insertion in the reference must be handled correctly.
        """
        template = (
            'CA-CGTG',
            ' AAC',
            ' ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'CAACGTG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testTwoInsertionsInReference(self):
        """
        Two insertions in the reference must be handled correctly.
        """
        template = (
            'CA-CGT-G',
            ' AAC',
            ' ???',
            '     TAG',
            '     ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'CAACGTAG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testTwoInsertionsAndOneDeletionInReference(self):
        """
        Two insertions and one deletion in the reference must be handled
        correctly.
        """
        template = (
            'CA-CGT-G',
            ' AA-G',
            ' ??-?',
            '     TAG',
            '     ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'CAAGTAG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testTwoInsertionsAndTwoDeletionsInReference(self):
        """
        Two insertions and two deletions in the reference must be handled
        correctly.
        """
        template = (
            'CA-CGT-G',
            ' AA--',
            ' ??--',
            '     TAG',
            '     ???',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'CAATAG',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.7,
                                 ignoreQuality=self.ignoreQuality))

    def testGeneiousExamplesNoTie(self):
        """
        Test the no-tied counts example from
        https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
        """
        template = (
            'ACGT',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  G',
            '  ?',
            '  G',
            '  ?',
            '  G',
            '  ?',
            '  T',
            '  ?',
        )

        with makeBAM(template) as (reference, bamFilename):
            for expected, threshold in ('A', 0.4), ('R', 0.7), ('D', 0.95):
                self.assertEqual(
                    f'AC{expected}T',
                    consensusFromBAM(bamFilename,
                                     reference=reference,
                                     threshold=threshold,
                                     ignoreQuality=self.ignoreQuality))

    def testGeneiousExamplesTie(self):
        """
        Test the tied counts example from
        https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
        """
        template = (
            'ACGT',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  A',
            '  ?',
            '  G',
            '  ?',
            '  G',
            '  ?',
            '  T',
            '  ?',
            '  T',
            '  ?',
        )

        with makeBAM(template) as (reference, bamFilename):
            for expected, threshold in ('A', 0.4), ('D', 0.7), ('D', 0.95):
                self.assertEqual(
                    f'AC{expected}T',
                    consensusFromBAM(bamFilename,
                                     reference=reference,
                                     threshold=threshold,
                                     ignoreQuality=self.ignoreQuality))

    def testOmicronEPE214Insertion(self):
        """
        Test that an amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) works as expected.

        The nucleotide sequence below can be obtained via:

        $ ncbi-fetch-id.py MN908947.3 > MN908947.3.fasta
        $ describe-genome.py --feature S --printNtSeq < MN908947.3.fasta | \
              filter-fasta.py --keepSites 630-673 --quiet | tail -n 1

        I then inserted the 9 nucleotide sequence GAGCCAGAA before the TGAT...
        starting at position 642.
        """
        template = (
            'TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
            '      AGTGCGGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTT',
            '      ??????????????????????????????????????????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'TAATTTAGTGCGGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 ignoreQuality=self.ignoreQuality))

    def testOmicronEPE214InsertionRightSideExact(self):
        """
        Test that an amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) works as expected when the insertion is the very
        last part of the read.
        """
        template = (
            'TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
            '      AGTGCGGAGCCAGAA',
            '      ???????????????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'TAATTTAGTGCGGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 ignoreQuality=self.ignoreQuality))

    def testOmicronEPE214InsertionLeftSideExact(self):
        """
        Test that an amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) works as expected when the insertion is the very
        beginning of the read.
        """
        template = (
            'TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
            '            GAGCCAGAATGATCT',
            '            ???????????????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'TAATTTAGTGCGGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 ignoreQuality=self.ignoreQuality))

    def testOmicronEPE214InsertionInTwoReads(self):
        """
        Test that an amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) works as expected, including when the insertion
        is inferred from two reads that each partially cover it. The
        nucleotide sequence below is obtained as in the test above.
        """
        template = (
            'TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
            '      AGTGCGGAGCCA',
            '      ????????????',
            '             AGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTT',
            '             ???????????????????????????????????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'TAATTTAGTGCGGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 ignoreQuality=self.ignoreQuality))

    def testOmicronEPE214PartialInsertion(self):
        """
        Test that a trailing part of the amino acid EPE sequence (here
        'AGCCAGAA') insertion into the SARS-CoV-2 spike nucleotide sequence
        that usually is found at location 642 (amino acid location 214) works
        as expected, including when the insertion is only partial (here we
        don't have the first nucleotide of the 9 nucleotide insertion).
        """
        template = (
            '--------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
            'AGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTT',
            '???????????????????????????????????',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'AGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 ignoreQuality=self.ignoreQuality))


@skipUnless(samtoolsInstalled(), 'samtools is not installed')
class TestIgnoreQuality(TestCase, _Mixin):
    """
    Test making majority consensuses with quality ignored.
    """
    ignoreQuality = True


@skipUnless(samtoolsInstalled(), 'samtools is not installed')
class TestWithQuality(TestCase, _Mixin):
    """
    Test making majority consensuses with quality considered.
    """
    ignoreQuality = False

    def testHighQualityDominates(self):
        """
        If one read has a very high quality base (here 'C', quality ']' = 60)
        that base should take precedence over two other reads that agree with
        each other but which have lower ('5' = 30) quality.
        """
        template = (
            'ACGT',
            '  A',
            '  5',
            '  A',
            '  5',
            '  C',
            '  ]',
        )

        with makeBAM(template) as (reference, bamFilename):
            self.assertEqual(
                'ACCT',
                consensusFromBAM(bamFilename,
                                 reference=reference,
                                 threshold=0.5,
                                 ignoreQuality=self.ignoreQuality))

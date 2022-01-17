from unittest import TestCase, skipUnless
from contextlib import contextmanager
from pathlib import Path
from tempfile import mkdtemp

from dark.consensus import consensusFromBAM
from dark.process import Executor
from dark.reads import DNARead
from dark.sam import samtoolsInstalled


@contextmanager
def makeBAM(template):
    """
    A context manager decorator to make a simple BAM file from a template.

    @param template: An iterable of C{str} sequences. The first will be treated
        as the reference, and then subsequent pairs (if any) will be treated as
        read and quality strings. Reads and quality strings can be indented
        with spaces to show where the read aligns with the reference (as would
        be produced by a tool like bowtie2 or bwa).
    @return: A context manager that produces a reference sequence and the name
        of the BAM file.
    """
    e = Executor()
    dirname = mkdtemp(prefix='test-consensus-')
    refId = 'ref-id'
    try:
        samFile = Path(dirname) / 'file.sam'
        bamFile = Path(dirname) / 'file.bam'
        reference = DNARead(refId, template[0])
        assert (len(template) % 2) == 1
        nSeqs = (len(template) - 1) >> 1

        with open(samFile, 'w') as fp:
            print(f'@SQ\tSN:{refId}\tLN:{len(template[0])}', file=fp)

            for count in range(nSeqs):
                read = template[count * 2 + 1].rstrip()
                quality = template[count * 2 + 2].rstrip()
                assert len(read) == len(quality)
                query = read.lstrip()
                quality = quality.lstrip()
                matchOffset = len(read) - len(query)
                print('\t'.join(map(str, (
                    f'read{count}',  # QUERY ID
                    0,  # FLAGS
                    refId,  # REF ID
                    matchOffset + 1,  # POS
                    30,  # MAPQ
                    f'{len(query)}M',  # CIGAR
                    '=',  # MRNM (Mate reference name)
                    1,  # MPOS (Mate position)
                    0,  # ISIZE
                    query,
                    quality,
                ))), file=fp)

        e.execute(f'samtools view -b -o {str(bamFile)!r} {str(samFile)!r}')
        e.execute(f'samtools index {str(bamFile)!r}')
        yield (reference, bamFile)
    finally:
        e.execute(f'rm -fr {dirname!r}')


class _Mixin:
    """
    Common (i.e., with and without considering quality) tests for consensuses
    making.
    """
    def testNoReadsReference(self):
        """
        If no reads are present and resolution of no-coverage bases is the
        reference sequence, the reference should be returned.
        """
        template = ('ACGTTCCG',)

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                template[0],
                consensusFromBAM(bamFilename, reference,
                                 noCoverage='reference',
                                 ignoreQuality=self.ignoreQuality))

    def testNoReadsN(self):
        """
        If no reads are present and resolution of no-coverage bases is 'N'
        (or '?', etc) a sequence of Ns (or ?s, etc) should be returned.
        """
        template = (
            'ACGTTCCG',
        )

        with makeBAM(template) as data:
            reference, bamFilename = data
            for char in 'N?':
                self.assertEqual(
                    char * len(template[0]),
                    consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                template[0],
                consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            for char in 'N?':
                self.assertEqual(
                    'AC' + char * 3 + 'CCG',
                    consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            for char in 'N?':
                self.assertEqual(
                    'XX' + char * 3 + 'XXX',
                    consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                template[0],
                consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                'ACAAACCG',
                consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                'ACAAA+CG',
                consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                'NNAAA+NN',
                consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                'ACAT',
                consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                'ACMT',
                consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            for expected, threshold in ('A', 0.4), ('R', 0.7), ('D', 0.95):
                self.assertEqual(
                    f'AC{expected}T',
                    consensusFromBAM(bamFilename, reference,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            for expected, threshold in ('A', 0.4), ('D', 0.7), ('D', 0.95):
                self.assertEqual(
                    f'AC{expected}T',
                    consensusFromBAM(bamFilename, reference,
                                     threshold=threshold,
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

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(
                'ACCT',
                consensusFromBAM(bamFilename, reference,
                                 threshold=0.5,
                                 ignoreQuality=self.ignoreQuality))

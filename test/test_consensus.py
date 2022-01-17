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

        # e.execute(f'cat {samFile}')
        e.execute(f'samtools view -b -o {str(bamFile)!r} {str(samFile)!r}')
        e.execute(f'samtools index {str(bamFile)!r}')
        yield (reference, bamFile)
    finally:
        e.execute(f'rm -fr {dirname!r}')


@skipUnless(samtoolsInstalled(), 'samtools is not installed')
class TestMajority(TestCase):
    """
    Test making majority consensuses.
    """
    def testNoReadsReference(self):
        """
        If no reads are present and resolution of no-coverage bases is the
        reference sequence, the reference should be returned.
        """
        template = ('ACGTTCCG',)

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(template[0],
                             consensusFromBAM(bamFilename, reference,
                                              strategy='majority',
                                              noCoverage='reference'))

    def testNoReadsN(self):
        """
        If no reads are present and resolution of no-coverage bases is 'N'
        (or '?', etc) a sequence of Ns (or ?s, etc) should be returned.
        """
        template = ('ACGTTCCG',)

        with makeBAM(template) as data:
            reference, bamFilename = data
            for char in 'N?':
                self.assertEqual(char * len(template[0]),
                                 consensusFromBAM(bamFilename, reference,
                                                  strategy='majority',
                                                  noCoverage=char))

    def testLowReadsReference(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        bases is the reference sequence, the reference should be returned.
        """
        template = ('ACGTTCCG',
                    '  GTT',
                    '  ???',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(template[0],
                             consensusFromBAM(bamFilename, reference,
                                              strategy='majority',
                                              lowCoverage='reference',
                                              minCoverage=2))

    def testLowReadsCharNoCoverageConsensus(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        bases is 'N' (or '? etc), then Ns (or ?s, etc) should be returned in
        the low coverage sites, and the reference in the sites with no
        coverage.
        """
        template = ('ACGTTCCG',
                    '  GTT',
                    '  ???',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            for char in 'N?':
                self.assertEqual('AC' + char * 3 + 'CCG',
                                 consensusFromBAM(bamFilename, reference,
                                                  strategy='majority',
                                                  lowCoverage=char,
                                                  minCoverage=2))

    def testLowReadsCharNoCoverageX(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        bases is 'N' (or '? etc), then Ns (or ?s, etc) should be returned in
        the low coverage sites, and (for example) 'X' in the sites with no
        coverage.
        """
        template = ('ACGTTCCG',
                    '  GTT',
                    '  ???',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            for char in 'N?':
                self.assertEqual('XX' + char * 3 + 'XXX',
                                 consensusFromBAM(bamFilename, reference,
                                                  strategy='majority',
                                                  noCoverage='X',
                                                  lowCoverage=char,
                                                  minCoverage=2))

    def testOneReadMatchingPartOfTheReference(self):
        """
        If one read is present and it matches part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        reference should be returned.
        """
        template = ('ACGTTCCG',
                    '  GTT',
                    '  ???',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual(template[0],
                             consensusFromBAM(bamFilename, reference,
                                              strategy='majority',
                                              noCoverage='reference'))

    def testOneReadDifferingFromPartOfTheReference(self):
        """
        If one read is present and it differs from part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        expected hybrid of the read and the reference should be returned.
        """
        template = ('ACGTTCCG',
                    '  AAA',
                    '  ???',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual('ACAAACCG',
                             consensusFromBAM(bamFilename, reference,
                                              strategy='majority',
                                              noCoverage='reference'))

    def testTwoReadsDifferingFromPartOfTheReferenceSomeLowCoverage(self):
        """
        If two reads are present and they differ from part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        expected hybrid of the read and the reference should be returned, with
        the part of the reads that has insufficient coverage returning the
        low-coverage symbol (here '+').
        """
        template = ('ACGTTCCG',
                    '  AAA',
                    '  ???',
                    '  AAAT',
                    '  ????',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual('ACAAA+CG',
                             consensusFromBAM(bamFilename, reference,
                                              strategy='majority',
                                              minCoverage=2,
                                              noCoverage='reference',
                                              lowCoverage='+'))

    def testTwoReadsDifferingFromPartOfTheReferenceLowAndNoCoverage(self):
        """
        If two reads are present and they differ from part of the reference and
        resolution of no-coverage bases is 'N', the expected hybrid of the read
        and the reference should be returned, with the part of the reads that
        has insufficient coverage returning the low-coverage symbol (here '+').
        """
        template = ('ACGTTCCG',
                    '  AAA',
                    '  ???',
                    '  AAAT',
                    '  ????',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual('NNAAA+NN',
                             consensusFromBAM(bamFilename, reference,
                                              strategy='majority',
                                              minCoverage=2,
                                              noCoverage='N',
                                              lowCoverage='+'))

    def testSimpleMajority(self):
        """
        If three reads result in a majority base at a site, that base should
        be in the consensus.
        """
        template = ('ACGT',
                    '  A',
                    '  ?',
                    '  A',
                    '  ?',
                    '  C',
                    '  ?',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual('ACAT',
                             consensusFromBAM(bamFilename, reference,
                                              threshold=0.5,
                                              strategy='majority'))

    def testSimpleMajorityBelowThreshold(self):
        """
        If conflicting reads at a site do not give a simple (above threshold)
        majority, the ambiguous code should be in the consensus.
        """
        template = ('ACGT',
                    '  A',
                    '  ?',
                    '  A',
                    '  ?',
                    '  C',
                    '  ?',
                    )

        with makeBAM(template) as data:
            reference, bamFilename = data
            self.assertEqual('ACMT',
                             consensusFromBAM(bamFilename, reference,
                                              threshold=0.7,
                                              strategy='majority'))

from contextlib import contextmanager
from pathlib import Path
from tempfile import mkdtemp

from dark.cigar import makeCigar
from dark.process import Executor
from dark.reads import DNARead
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

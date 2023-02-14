from contextlib import contextmanager
from pathlib import Path
from tempfile import mkdtemp

from dark.cigar import makeCigar
from dark.process import Executor
from dark.reads import DNARead, Reads
from dark.utils import matchOffset

# From https://samtools.github.io/hts-specs/SAMv1.pdf
CINS, CDEL, CMATCH = tuple("IDM")

REF_ID = "ref-id"


@contextmanager
def makeBAM(template, bamReferences=None, fastaReferences=None):
    """
    A context manager decorator to make a simple BAM file from a template.
    Note that this code invokes samtools.

    @param template: An iterable of C{str} sequences. The first will be treated
        as the reference, and then subsequent pairs (if any) will be treated as
        read and quality strings. Reads and quality strings can be indented
        with spaces to show where the read aligns with the reference.
    @return: A context manager that produces a 2-tuple containing the reference
        C{DNARead} instance and the C{Path} of the BAM file.
    """
    if len(template) % 2 != 1:
        raise ValueError(
            "The template must have an odd number of strings, specifying the "
            "reference sequence, then zero or more read/quality pairs."
        )

    leftPaddedReference = template[0]
    templateSequence = leftPaddedReference.lstrip().replace("-", "")

    if bamReferences is None:
        matchedReference = DNARead(REF_ID, templateSequence)
        bamReferences = Reads([matchedReference])
    else:
        matchedReference = bamReferences[0]
        # Sanity check: The first BAM reference must have the same sequence
        # as the template.
        assert matchedReference.sequence == templateSequence
        bamReferences = Reads(bamReferences)

    fastaReferences = Reads(
        bamReferences if fastaReferences is None else fastaReferences
    )

    nSeqs = (len(template) - 1) >> 1
    dirname = mkdtemp(prefix="test-consensus-")
    e = Executor()

    try:
        fastaFile = Path(dirname) / "references.fasta"
        samFile = Path(dirname) / "file.sam"
        bamFile = Path(dirname) / "file.bam"

        fastaReferences.save(fastaFile)

        with open(samFile, "w") as fp:
            for reference in bamReferences:
                print(f"@SQ\tSN:{reference.id}\tLN:{len(reference)}", file=fp)

            for count in range(nSeqs):
                leftPaddedQuery = template[count * 2 + 1].rstrip()
                leftPaddedQuality = template[count * 2 + 2].rstrip()
                assert len(leftPaddedQuery) == len(leftPaddedQuality)
                query = leftPaddedQuery.lstrip()
                quality = leftPaddedQuality.lstrip()
                queryNoGaps = qualityNoGaps = ""
                for queryBase, qualityBase in zip(query, quality):
                    if queryBase != "-":
                        queryNoGaps += queryBase
                        qualityNoGaps += qualityBase

                print(
                    "\t".join(
                        map(
                            str,
                            (
                                f"read{count}",  # QNAME (query name)
                                0,  # FLAGS
                                matchedReference.id,  # RNAME (reference name)
                                matchOffset(leftPaddedReference, leftPaddedQuery) + 1,
                                30,  # MAPQ (mapping quality)
                                makeCigar(
                                    leftPaddedReference, leftPaddedQuery
                                ),  # CIGAR
                                "*",  # MRNM (mate reference name)
                                0,  # MPOS (mate position)
                                0,  # ISIZE (insert size)
                                queryNoGaps,  # SEQ
                                qualityNoGaps,  # QUAL
                            ),
                        )
                    ),
                    file=fp,
                )

        e.execute(
            f"samtools sort -O BAM --write-index -o {str(bamFile)!r} "
            f"{str(samFile)!r}"
        )
        yield (fastaFile, bamFile)
    finally:
        # import sys; print(f'{samFile}', file=sys.stderr)
        e.execute(f"rm -fr {dirname!r}")

from contextlib import contextmanager
from os import close, unlink, write
from tempfile import mkstemp
from unittest import TestCase

from dark.sam import UnknownReference, UnspecifiedReference, coverageDepth


@contextmanager
def dataFile(data):
    """
    Create a temporary file containing C{data} and yield its filename,
    removing it on exit.
    """
    fd, filename = mkstemp()
    write(fd, data.encode("utf-8"))
    close(fd)
    yield filename
    unlink(filename)


def samData(*read_lines, references=(("ref1", 100),)):
    """
    Build a SAM string from header references and alignment lines.

    @param read_lines: Alignment lines (space-separated, will be tab-converted).
    @param references: An iterable of (name, length) pairs for @SQ headers.
    @return: A C{str} SAM-formatted string.
    """
    lines = [f"@SQ SN:{name} LN:{length}" for name, length in references]
    lines.extend(read_lines)
    return "\n".join(lines).replace(" ", "\t")


class TestCoverageDepth(TestCase):
    """
    Tests for dark.sam.coverageDepth.
    """

    # --- reference explicitly given ---

    def testUnknownReference(self):
        """
        Passing an unknown reference name must raise UnknownReference.
        """
        data = samData(references=[("ref1", 10)])
        with dataFile(data) as filename:
            self.assertRaises(UnknownReference, coverageDepth, filename, "unknown")

    def testNoReads(self):
        """
        A SAM file with no aligned reads must return all-zero depth.
        """
        data = samData(references=[("ref1", 10)])
        with dataFile(data) as filename:
            ref, result = coverageDepth(filename, "ref1")
        self.assertEqual("ref1", ref)
        self.assertEqual([0] * 10, result)

    def testSingleReadFullyCovering(self):
        """
        A single read covering the whole reference must give depth 1 everywhere.
        """
        # Read at 1-based position 1, CIGAR 10M → covers 0-based positions 0-9.
        data = samData(
            "read1 0 ref1 1 60 10M * 0 0 ACGTACGTAC *",
            references=[("ref1", 10)],
        )
        with dataFile(data) as filename:
            ref, result = coverageDepth(filename, "ref1")
        self.assertEqual("ref1", ref)
        self.assertEqual([1] * 10, result)

    def testSingleReadPartialCoverage(self):
        """
        A read starting at position 3 (1-based) covering 4 bases must give
        depth 1 at 0-based positions 2-5 and 0 elsewhere.
        """
        data = samData(
            "read1 0 ref1 3 60 4M * 0 0 ACGT *",
            references=[("ref1", 10)],
        )
        with dataFile(data) as filename:
            _, result = coverageDepth(filename, "ref1")
        self.assertEqual([0, 0, 1, 1, 1, 1, 0, 0, 0, 0], result)

    def testTwoOverlappingReads(self):
        """
        Two overlapping reads must produce depth 2 in the overlap region.
        """
        # read1: positions 1-5 (0-based 0-4), read2: positions 3-7 (0-based 2-6).
        data = samData(
            "read1 0 ref1 1 60 5M * 0 0 ACGTA *",
            "read2 0 ref1 3 60 5M * 0 0 GTACG *",
            references=[("ref1", 10)],
        )
        with dataFile(data) as filename:
            _, result = coverageDepth(filename, "ref1")
        self.assertEqual([1, 1, 2, 2, 2, 1, 1, 0, 0, 0], result)

    def testUnmappedReadsIgnored(self):
        """
        Unmapped reads (FLAG 4) must not contribute to coverage depth.
        """
        data = samData(
            "read1 4 * 0 0 * * 0 0 ACGT *",
            references=[("ref1", 10)],
        )
        with dataFile(data) as filename:
            _, result = coverageDepth(filename, "ref1")
        self.assertEqual([0] * 10, result)

    def testMultipleReferences(self):
        """
        Only reads mapping to the requested reference must be counted.
        """
        data = samData(
            "read1 0 ref1 1 60 5M * 0 0 ACGTA *",
            "read2 0 ref2 1 60 5M * 0 0 ACGTA *",
            references=[("ref1", 10), ("ref2", 10)],
        )
        with dataFile(data) as filename:
            _, depth1 = coverageDepth(filename, "ref1")
            _, depth2 = coverageDepth(filename, "ref2")
        self.assertEqual([1, 1, 1, 1, 1, 0, 0, 0, 0, 0], depth1)
        self.assertEqual([1, 1, 1, 1, 1, 0, 0, 0, 0, 0], depth2)

    def testDeletionInCigar(self):
        """
        Deleted reference bases (CIGAR D operation) must not inflate depth.
        """
        # 2M1D2M at pos 1: covers 0-based positions 0,1 (match), 2 (deleted,
        # not counted), 3,4 (match). Depth should be 1 at 0,1,3,4 and 0 at 2.
        data = samData(
            "read1 0 ref1 1 60 2M1D2M * 0 0 ACGT *",
            references=[("ref1", 10)],
        )
        with dataFile(data) as filename:
            _, result = coverageDepth(filename, "ref1")
        self.assertEqual([1, 1, 0, 1, 1, 0, 0, 0, 0, 0], result)

    def testReturnLengthMatchesReference(self):
        """
        The returned list length must equal the reference sequence length.
        """
        data = samData(references=[("ref1", 50)])
        with dataFile(data) as filename:
            _, result = coverageDepth(filename, "ref1")
        self.assertEqual(50, len(result))

    # --- reference auto-detected (not given) ---

    def testAutoDetectSingleReference(self):
        """
        When no reference is given and all reads map to one reference,
        that reference must be auto-detected and returned.
        """
        data = samData(
            "read1 0 ref1 1 60 5M * 0 0 ACGTA *",
            "read2 0 ref1 3 60 5M * 0 0 GTACG *",
            references=[("ref1", 10)],
        )
        with dataFile(data) as filename:
            ref, result = coverageDepth(filename)
        self.assertEqual("ref1", ref)
        self.assertEqual([1, 1, 2, 2, 2, 1, 1, 0, 0, 0], result)

    def testAutoDetectFailsOnMultipleReferences(self):
        """
        When no reference is given and reads map to multiple references,
        UnspecifiedReference must be raised.
        """
        data = samData(
            "read1 0 ref1 1 60 5M * 0 0 ACGTA *",
            "read2 0 ref2 1 60 5M * 0 0 ACGTA *",
            references=[("ref1", 10), ("ref2", 10)],
        )
        with dataFile(data) as filename:
            self.assertRaises(UnspecifiedReference, coverageDepth, filename)

    def testAutoDetectFailsOnNoMappedReads(self):
        """
        When no reference is given and there are no mapped reads,
        UnspecifiedReference must be raised.
        """
        data = samData(
            "read1 4 * 0 0 * * 0 0 ACGT *",
            references=[("ref1", 10)],
        )
        with dataFile(data) as filename:
            self.assertRaises(UnspecifiedReference, coverageDepth, filename)

    def testAutoDetectIgnoresHeaderOnlyReferences(self):
        """
        When no reference is given, references listed only in the SAM header
        but with no aligned reads must not trigger a multiple-reference error.
        """
        # Header mentions ref1 and ref2, but only ref1 has aligned reads.
        data = samData(
            "read1 0 ref1 1 60 5M * 0 0 ACGTA *",
            references=[("ref1", 10), ("ref2", 10)],
        )
        with dataFile(data) as filename:
            ref, result = coverageDepth(filename)
        self.assertEqual("ref1", ref)
        self.assertEqual([1, 1, 1, 1, 1, 0, 0, 0, 0, 0], result)

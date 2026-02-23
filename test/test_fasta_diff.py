import importlib
import sys
from io import StringIO
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from dark.reads import DNARead

# The script filename contains a hyphen so we load it via importlib.
# Temporarily disable bytecode caching to avoid creating __pycache__ in bin/
# which breaks setuptools data-files glob.
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
with patch.object(sys, "dont_write_bytecode", True):
    fasta_diff = importlib.import_module("fasta-diff")


C = fasta_diff.C
Difference = fasta_diff.Difference
findDifferences = fasta_diff.findDifferences
formatBasePair = fasta_diff.formatBasePair
pickSeq = fasta_diff.pickSeq
printAlignment = fasta_diff.printAlignment
printDiffTable = fasta_diff.printDiffTable
printSummary = fasta_diff.printSummary
main = fasta_diff.main


class TestPickSeq(TestCase):
    """
    Test the pickSeq function.
    """

    def testFirstSequenceByDefault(self):
        """
        When no name is given, the first sequence must be returned.
        """
        reads = [DNARead("seq1", "ACGT"), DNARead("seq2", "TGCA")]
        result = pickSeq(reads, None, "test.fa")
        self.assertEqual("seq1", result.id)

    def testPickByName(self):
        """
        When a name is given, the matching sequence must be returned.
        """
        reads = [DNARead("seq1", "ACGT"), DNARead("seq2", "TGCA")]
        result = pickSeq(reads, "seq2", "test.fa")
        self.assertEqual("seq2", result.id)

    def testMissingNameExits(self):
        """
        When a name is given that does not exist, SystemExit must be raised.
        """
        reads = [DNARead("seq1", "ACGT")]
        self.assertRaises(SystemExit, pickSeq, reads, "missing", "test.fa")

    def testEmptySeqsExits(self):
        """
        When no sequences are available, SystemExit must be raised.
        """
        self.assertRaises(SystemExit, pickSeq, [], None, "test.fa")


class TestFindDifferences(TestCase):
    """
    Test the findDifferences function.
    """

    def testIdenticalSequences(self):
        """
        Two identical sequences must produce no differences.
        """
        self.assertEqual([], findDifferences("ACGT", "ACGT"))

    def testSingleMismatch(self):
        """
        A single mismatch must be found.
        """
        diffs = findDifferences("ACGT", "ACAT")
        self.assertEqual(1, len(diffs))
        self.assertEqual(
            Difference(
                pos=2, refPos=2, qryPos=2, refBase="G", qryBase="A", kind="mismatch"
            ),
            diffs[0],
        )

    def testMultipleMismatches(self):
        """
        Multiple mismatches must all be found.
        """
        diffs = findDifferences("AAAA", "TTTT")
        self.assertEqual(4, len(diffs))
        for d in diffs:
            self.assertEqual("mismatch", d.kind)

    def testDeletion(self):
        """
        A gap in the second sequence must be reported as a deletion.
        """
        diffs = findDifferences("ACGT", "AC-T")
        self.assertEqual(1, len(diffs))
        d = diffs[0]
        self.assertEqual("deletion", d.kind)
        self.assertEqual(2, d.refPos)
        self.assertIsNone(d.qryPos)
        self.assertEqual("G", d.refBase)
        self.assertEqual("-", d.qryBase)

    def testInsertion(self):
        """
        A gap in the first sequence must be reported as an insertion.
        """
        diffs = findDifferences("AC-T", "ACGT")
        self.assertEqual(1, len(diffs))
        d = diffs[0]
        self.assertEqual("insertion", d.kind)
        self.assertIsNone(d.refPos)
        self.assertEqual(2, d.qryPos)
        self.assertEqual("-", d.refBase)
        self.assertEqual("G", d.qryBase)

    def testEmptySequences(self):
        """
        Two empty sequences must produce no differences.
        """
        self.assertEqual([], findDifferences("", ""))

    def testRefPositionSkipsGaps(self):
        """
        The refPos must not count gap characters in the first sequence.
        """
        # aln: A - C G
        # ref: 0   1 2  (gaps skipped)
        diffs = findDifferences("A-CG", "AACG")
        self.assertEqual(1, len(diffs))
        self.assertIsNone(diffs[0].refPos)
        self.assertEqual(1, diffs[0].qryPos)

    def testQryPositionSkipsGaps(self):
        """
        The qryPos must not count gap characters in the second sequence.
        """
        diffs = findDifferences("AACG", "A-CG")
        self.assertEqual(1, len(diffs))
        self.assertEqual(1, diffs[0].refPos)
        self.assertIsNone(diffs[0].qryPos)

    def testMixedDifferences(self):
        """
        A mix of mismatches, insertions, and deletions must all be found.
        """
        #                          D       I  M
        diffs = findDifferences("ATC-GA", "A-CGTA")
        kinds = [d.kind for d in diffs]
        self.assertIn("mismatch", kinds)
        self.assertIn("insertion", kinds)
        self.assertIn("deletion", kinds)


class TestFormatBasePair(TestCase):
    """
    Test the formatBasePair function.
    """

    def setUp(self):
        C.disable()

    def testMatchingBases(self):
        """
        Matching bases must produce a space as the diff indicator.
        """
        _, _, indicator = formatBasePair("A", "A")
        self.assertEqual(" ", indicator)

    def testMatchingBasesReturnBases(self):
        """
        Matching bases must return the base character in both sequence slots.
        """
        fa, fb, _ = formatBasePair("A", "A")
        self.assertIn("A", fa)
        self.assertIn("A", fb)

    def testMismatch(self):
        """
        Mismatched bases must produce a '*' diff indicator.
        """
        _, _, indicator = formatBasePair("A", "T")
        self.assertEqual("*", indicator)

    def testMismatchReturnsBothBases(self):
        """
        Mismatched bases must return each base in its respective sequence slot.
        """
        fa, fb, _ = formatBasePair("A", "T")
        self.assertIn("A", fa)
        self.assertIn("T", fb)

    def testInsertionIndicator(self):
        """
        A gap in the first sequence must produce a '+' diff indicator.
        """
        _, _, indicator = formatBasePair("-", "A")
        self.assertEqual("+", indicator)

    def testDeletionIndicator(self):
        """
        A gap in the second sequence must produce a '-' diff indicator.
        """
        _, _, indicator = formatBasePair("A", "-")
        self.assertEqual("-", indicator)


class TestPrintSummary(TestCase):
    """
    Test the printSummary function.
    """

    def setUp(self):
        C.disable()

    def testNoDifferences(self):
        """
        When there are no differences, 100% identity must be reported.
        """
        captured = StringIO()
        sys.stdout = captured
        try:
            printSummary([], 10, 10, 10)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("100.00%", output)
        self.assertIn("Total diffs:        0", output)

    def testWithDifferences(self):
        """
        When there are differences, the correct counts must be reported.
        """
        diffs = [
            Difference(0, 0, 0, "A", "T", "mismatch"),
            Difference(1, 1, None, "C", "-", "deletion"),
            Difference(2, None, 1, "-", "G", "insertion"),
        ]
        captured = StringIO()
        sys.stdout = captured
        try:
            printSummary(diffs, 10, 10, 10)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("Mismatches:       1", output)
        self.assertIn("Deletions:        1", output)
        self.assertIn("Insertions:       1", output)
        self.assertIn("Total diffs:        3", output)

    def testAlignmentLengthShown(self):
        """
        When the alignment length differs from the sequence lengths, it must
        be printed.
        """
        captured = StringIO()
        sys.stdout = captured
        try:
            printSummary([], 8, 10, 12)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("Alignment len:", output)

    def testAlignmentLengthNotShownWhenEqual(self):
        """
        When the alignment length equals both sequence lengths, the alignment
        length line must not be printed.
        """
        captured = StringIO()
        sys.stdout = captured
        try:
            printSummary([], 10, 10, 10)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertNotIn("Alignment len:", output)

    def testZeroAlignmentLength(self):
        """
        When the alignment length is zero, identity must default to 100%.
        """
        captured = StringIO()
        sys.stdout = captured
        try:
            printSummary([], 0, 0, 0)
        finally:
            sys.stdout = sys.__stdout__
        self.assertIn("100.00%", captured.getvalue())


class TestPrintDiffTable(TestCase):
    """
    Test the printDiffTable function.
    """

    def setUp(self):
        C.disable()

    def testNoDifferences(self):
        """
        When there are no differences, nothing must be printed.
        """
        captured = StringIO()
        sys.stdout = captured
        try:
            printDiffTable([])
        finally:
            sys.stdout = sys.__stdout__
        self.assertEqual("", captured.getvalue())

    def testSingleDifference(self):
        """
        A single difference must appear in the table.
        """
        diffs = [Difference(5, 5, 5, "A", "T", "mismatch")]
        captured = StringIO()
        sys.stdout = captured
        try:
            printDiffTable(diffs)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("mismatch", output)
        self.assertIn("Differences", output)

    def testCompactLimit(self):
        """
        When the number of differences exceeds the compact limit, a truncation
        message must be shown.
        """
        diffs = [Difference(i, i, i, "A", "T", "mismatch") for i in range(10)]
        captured = StringIO()
        sys.stdout = captured
        try:
            printDiffTable(diffs, compactLimit=3)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("showing first 3", output)
        self.assertIn("... and 7 more", output)

    def testInsertionRowShowsDash(self):
        """
        An insertion (None refPos) must display an em-dash in the Seq1 column.
        """
        diffs = [Difference(5, None, 4, "-", "G", "insertion")]
        captured = StringIO()
        sys.stdout = captured
        try:
            printDiffTable(diffs)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("—", output)
        self.assertIn("insertion", output)

    def testDeletionRowShowsDash(self):
        """
        A deletion (None qryPos) must display an em-dash in the Seq2 column.
        """
        diffs = [Difference(5, 4, None, "G", "-", "deletion")]
        captured = StringIO()
        sys.stdout = captured
        try:
            printDiffTable(diffs)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("—", output)
        self.assertIn("deletion", output)


class TestPrintAlignment(TestCase):
    """
    Test the printAlignment function.
    """

    def setUp(self):
        C.disable()

    def testIdenticalSequences(self):
        """
        When there are no differences, a message saying they are identical
        must be printed.
        """
        captured = StringIO()
        sys.stdout = captured
        try:
            printAlignment("ACGT", "ACGT", "s1", "s2", [])
        finally:
            sys.stdout = sys.__stdout__
        self.assertIn("Sequences are identical", captured.getvalue())

    def testDifferencesShown(self):
        """
        Differences must appear in the alignment output.
        """
        diffs = [Difference(2, 2, 2, "G", "A", "mismatch")]
        captured = StringIO()
        sys.stdout = captured
        try:
            printAlignment("ACGT", "ACAT", "s1", "s2", diffs, context=2)
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()
        self.assertIn("s1", output)
        self.assertIn("s2", output)

    def testIdentityMessageShowsLength(self):
        """
        The "Sequences are identical" message must include the alignment length.
        """
        captured = StringIO()
        sys.stdout = captured
        try:
            printAlignment("ACGTACGT", "ACGTACGT", "s1", "s2", [])
        finally:
            sys.stdout = sys.__stdout__
        self.assertIn("8 bp", captured.getvalue())

    def testDistantDiffsShowSeparator(self):
        """
        When two differences are far enough apart that their context windows
        do not overlap, a '...' separator must appear between the regions.
        """
        # Position the differences 40 apart so context=2 windows never merge.
        aln1 = "T" + "A" * 39 + "T"
        aln2 = "A" * 40 + "A"
        diffs = [
            Difference(0, 0, 0, "T", "A", "mismatch"),
            Difference(40, 40, 40, "T", "A", "mismatch"),
        ]
        captured = StringIO()
        sys.stdout = captured
        try:
            printAlignment(aln1, aln2, "s1", "s2", diffs, context=2)
        finally:
            sys.stdout = sys.__stdout__
        self.assertIn("...", captured.getvalue())


class TestMain(TestCase):
    """
    Test the main function via command-line argument simulation.
    """

    def _runMain(self, argv, reads1, reads2):
        """
        Run main() with mocked FastaReads and captured stdout.

        @param argv: A C{list} of command-line arguments (including the
            program name as the first element).
        @param reads1: A C{list} of C{Read} instances for file1.
        @param reads2: A C{list} of C{Read} instances for file2.
        @return: The C{str} captured stdout.
        """
        captured = StringIO()
        with patch.object(fasta_diff, "FastaReads") as mock_fr:
            mock_fr.side_effect = [reads1, reads2]
            sys.stdout = captured
            sys.argv = argv
            try:
                main()
            finally:
                sys.stdout = sys.__stdout__
        return captured.getvalue()

    def testIdenticalSequences(self):
        """
        Running main with two identical sequences must report 100% identity.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor"],
            [DNARead("s1", "ACGTACGT")],
            [DNARead("s2", "ACGTACGT")],
        )
        self.assertIn("100.00%", output)
        self.assertIn("Sequences are identical", output)

    def testSingleMismatch(self):
        """
        Running main with a single mismatch must report it.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACAT")],
        )
        self.assertIn("Mismatches:       1", output)
        self.assertIn("mismatch", output)

    def testDifferentLengthsWithEdlib(self):
        """
        Running main with sequences of different lengths must trigger
        alignment.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--aligner", "edlib"],
            [DNARead("s1", "ACGTACGT")],
            [DNARead("s2", "ACGT")],
        )
        self.assertIn("aligning with edlib", output)
        self.assertIn("Summary", output)

    def testPickByName(self):
        """
        The --name1 and --name2 arguments must select the correct sequences.
        """
        output = self._runMain(
            [
                "fasta-diff",
                "a.fa",
                "b.fa",
                "--noColor",
                "--name1",
                "second",
                "--name2",
                "beta",
            ],
            [DNARead("first", "AAAA"), DNARead("second", "ACGT")],
            [DNARead("alpha", "TTTT"), DNARead("beta", "ACGT")],
        )
        self.assertIn("second", output)
        self.assertIn("beta", output)
        self.assertIn("100.00%", output)

    def testNoTable(self):
        """
        The --noTable flag must suppress the diff table output.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--noTable"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACAT")],
        )
        self.assertIn("Summary", output)
        self.assertNotIn("mismatch", output)

    def testMaxDiffsExceeded(self):
        """
        When --maxDiffs is exceeded, only the summary must be printed, without
        the diff table or alignment.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--maxDiffs", "2"],
            [DNARead("s1", "AAAA")],
            [DNARead("s2", "TTTT")],
        )
        self.assertIn("Summary", output)
        self.assertIn("Total diffs:        4", output)
        # The diff table and alignment should be suppressed.
        self.assertNotIn("mismatch", output)
        self.assertNotIn("*", output)

    def testMaxDiffsNotExceeded(self):
        """
        When --maxDiffs is not exceeded, the diff table and alignment must
        be printed.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--maxDiffs", "5"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACAT")],
        )
        self.assertIn("mismatch", output)

    def testMaxDiffsMinusOne(self):
        """
        When --maxDiffs is -1 (the default), differences must always be
        printed regardless of count.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--maxDiffs", "-1"],
            [DNARead("s1", "AAAA")],
            [DNARead("s2", "TTTT")],
        )
        self.assertIn("mismatch", output)

    def testMaxDiffsZero(self):
        """
        When --maxDiffs is 0, any differences at all must suppress the diff
        table and alignment.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--maxDiffs", "0"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACAT")],
        )
        self.assertIn("Summary", output)
        self.assertNotIn("mismatch", output)

    def testMaxDiffsZeroIdentical(self):
        """
        When --maxDiffs is 0 and the sequences are identical, the identical
        message must still be printed (0 diffs <= 0 maxDiffs).
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--maxDiffs", "0"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACGT")],
        )
        self.assertIn("Sequences are identical", output)

    def testFullAlignment(self):
        """
        The --full flag must show the entire alignment without context
        windowing.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--full"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACAT")],
        )
        self.assertIn("s1", output)
        self.assertIn("s2", output)

    def testContextFlag(self):
        """
        The --context flag must be accepted and produce alignment output.
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--context", "1"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACAT")],
        )
        self.assertIn("s1", output)
        self.assertIn("mismatch", output)

    def testName1NotFound(self):
        """
        When --name1 specifies a name absent from file1, the script must exit
        with an error.
        """
        with patch.object(fasta_diff, "FastaReads") as mock_fr:
            mock_fr.side_effect = [
                [DNARead("s1", "ACGT")],
                [DNARead("s2", "ACGT")],
            ]
            sys.argv = [
                "fasta-diff",
                "a.fa",
                "b.fa",
                "--noColor",
                "--name1",
                "missing",
            ]
            self.assertRaises(SystemExit, main)

    def testMaxDiffsEqualToDiffs(self):
        """
        When --maxDiffs equals the number of differences exactly, the diff
        table and alignment must still be printed (boundary: <= not <).
        """
        output = self._runMain(
            ["fasta-diff", "a.fa", "b.fa", "--noColor", "--maxDiffs", "1"],
            [DNARead("s1", "ACGT")],
            [DNARead("s2", "ACAT")],
        )
        self.assertIn("mismatch", output)

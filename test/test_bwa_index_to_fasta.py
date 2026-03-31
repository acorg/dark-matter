import importlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest import TestCase, skipUnless
from unittest.mock import mock_open, patch

# The script filename contains a hyphen so we load it via importlib.
# Temporarily disable bytecode caching to avoid creating __pycache__ in bin/
# which breaks setuptools data-files glob.
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
with patch.object(sys, "dont_write_bytecode", True):
    bwa_to_fasta = importlib.import_module("bwa-index-to-fasta")

readAnn = bwa_to_fasta.readAnn
readAmb = bwa_to_fasta.readAmb
readPac = bwa_to_fasta.readPac
applyAmbiguities = bwa_to_fasta.applyAmbiguities
reconstructFasta = bwa_to_fasta.reconstructFasta

BWA_AVAILABLE = shutil.which("bwa") is not None


def basesToPac(*bases: str):
    """Pack a sequence of bases (strings) into BWA .pac bytes.

    Returns bytes suitable for passing to readPac via mock_open:
    the data bytes followed by the remainder byte.
    See for the .pac format see https://deepwiki.com/j-levy/bwa-gasal2/6.1-index-file-structure#packed-sequence-file-format-pac
    """
    # 0=A 1=C 2=G 3=T
    enc = {"A": 0, "C": 1, "G": 2, "T": 3}
    values = [enc[b] for b in bases]
    data: list[int] = []
    for i in range(0, len(values), 4):
        chunk = values[i : i + 4]
        byte = 0
        for j, v in enumerate(chunk):
            byte |= v << (6 - 2 * j)
        data.append(byte)
    remainder = len(values) % 4  # 0 means the last byte is full
    return bytes(data + [remainder])


class TestReadAnn(TestCase):
    """Tests for readAnn()."""

    def testSingleSequenceNoAnnotation(self):
        """A .ann file with one sequence and no annotation is parsed correctly."""
        name = "chr1"
        length = 50
        ann = f"{length} 1 0\n42 {name}\n0 {length} 0\n"
        with patch("builtins.open", mock_open(read_data=ann)):
            total, seqs = readAnn("fake.ann")

        self.assertEqual(length, total)
        self.assertEqual(1, len(seqs))
        self.assertEqual(name, seqs[0]["name"])
        self.assertEqual("", seqs[0]["anno"])
        self.assertEqual(0, seqs[0]["offset"])
        self.assertEqual(length, seqs[0]["length"])

    def testSingleSequenceWithAnnotation(self):
        """
        A .ann file where the name line has an annotation field is parsed correctly.
        """
        name = "seq1"
        annotation = "some annotation text"
        length = 100
        ann = f"{length} 1 0\n1 {name} {annotation}\n0 {length} 0\n"
        with patch("builtins.open", mock_open(read_data=ann)):
            total, seqs = readAnn("fake.ann")

        self.assertEqual(length, total)
        self.assertEqual(name, seqs[0]["name"])
        self.assertEqual(annotation, seqs[0]["anno"])
        self.assertEqual(0, seqs[0]["offset"])
        self.assertEqual(length, seqs[0]["length"])

    def testMultipleSequences(self):
        """A .ann file with two sequences returns both in order with correct offsets."""
        name1, len1 = "chrA", 80
        name2, anno2, len2 = "chrB", "chromosome B", 120
        totalLen = len1 + len2
        ann = (
            f"{totalLen} 2 0\n1 {name1}\n0 {len1} 0\n2 {name2} {anno2}\n{len1} "
            f"{len2} 3\n"
        )
        with patch("builtins.open", mock_open(read_data=ann)):
            total, seqs = readAnn("fake.ann")

        self.assertEqual(totalLen, total)
        self.assertEqual(2, len(seqs))

        self.assertEqual(name1, seqs[0]["name"])
        self.assertEqual("", seqs[0]["anno"])
        self.assertEqual(0, seqs[0]["offset"])
        self.assertEqual(len1, seqs[0]["length"])

        self.assertEqual(name2, seqs[1]["name"])
        self.assertEqual(anno2, seqs[1]["anno"])
        self.assertEqual(len1, seqs[1]["offset"])
        self.assertEqual(len2, seqs[1]["length"])

    def testTotalLengthReturnedFromHeader(self):
        """The first integer on the header line is returned as totalLength."""
        length = 999
        ann = f"{length} 1 0\n0 x\n0 {length} 0\n"
        with patch("builtins.open", mock_open(read_data=ann)):
            total, _ = readAnn("fake.ann")
        self.assertEqual(length, total)

    def testNullAnnotationTreatedAsEmpty(self):
        """BWA writes '(null)' when there's no annotation; this should be empty."""
        name = "chr1"
        length = 50
        ann = f"{length} 1 0\n0 {name} (null)\n0 {length} 0\n"
        with patch("builtins.open", mock_open(read_data=ann)):
            _, seqs = readAnn("fake.ann")
        self.assertEqual("", seqs[0]["anno"])


class TestReadAmb(TestCase):
    """Tests for readAmb()."""

    def testNoHoles(self):
        """A .amb file with zero holes returns an empty list."""
        expected = []
        amb = "100 1 0\n"
        with patch("builtins.open", mock_open(read_data=amb)):
            holes = readAmb("fake.amb")
        self.assertEqual(expected, holes)

    def testSingleHole(self):
        """A .amb file with one N-run is parsed correctly."""
        expected = [(10, 5, "N")]
        amb = "100 1 1\n10 5 N\n"
        with patch("builtins.open", mock_open(read_data=amb)):
            holes = readAmb("fake.amb")
        self.assertEqual(expected, holes)

    def testMultipleHoles(self):
        """A .amb file with multiple N-runs returns all of them in order."""
        expected = [(0, 4, "N"), (100, 10, "N"), (300, 2, "N")]
        amb = "500 2 3\n0 4 N\n100 10 N\n300 2 N\n"
        with patch("builtins.open", mock_open(read_data=amb)):
            holes = readAmb("fake.amb")
        self.assertEqual(expected, holes)

    def testHoleCharacterIsPreserved(self):
        """The ambiguity character recorded in .amb is stored as-is."""
        char = "X"
        amb = f"100 1 1\n5 3 {char}\n"
        with patch("builtins.open", mock_open(read_data=amb)):
            holes = readAmb("fake.amb")
        self.assertEqual(char, holes[0][2])

    def testMultipleHolesCharacterIsPreserved(self):
        """The ambiguity character recorded in .amb is stored as-is."""
        char1, char2 = "X", "N"
        amb = f"100 1 2\n5 3 {char1}\n11 4 {char2}\n"
        with patch("builtins.open", mock_open(read_data=amb)):
            holes = readAmb("fake.amb")
        self.assertEqual(char1, holes[0][2])
        self.assertEqual(char2, holes[1][2])


class TestReadPac(TestCase):
    """Tests for readPac()."""

    def testAllFourBaseTypes(self):
        """All four base types (A=00, C=01, G=10, T=11) decode correctly."""
        bases = list("ACGT")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)

    def testAllZeroBits(self):
        """A sequence of all A bases (zero bits) decodes correctly."""
        bases = list("AAAA")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)

    def testAllOneBits(self):
        """A sequence of all T bases (all one-bits) decodes correctly."""
        bases = list("TTTT")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)

    def testRemainderOne(self):
        """A single base (remainder=1) decodes correctly."""
        bases = list("T")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)

    def testRemainderTwo(self):
        """Two bases (remainder=2) decode correctly."""
        bases = list("CG")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)

    def testRemainderThree(self):
        """Three bases (remainder=3) decode correctly."""
        bases = list("GCA")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)

    def testMultipleFullBytes(self):
        """A sequence spanning two full data bytes decodes correctly."""
        bases = list("ACGTTGCA")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)

    def testMultipleBytesWithPartial(self):
        """Five bases (one full byte + one partial) decode correctly."""
        bases = list("ACGTA")
        pac = basesToPac(*bases)
        with patch("builtins.open", mock_open(read_data=pac)):
            seq = readPac("fake.pac")
        self.assertEqual(bases, seq)


class TestApplyAmbiguities(TestCase):
    """Tests for applyAmbiguities()."""

    def testNoHoles(self):
        """An empty holes list leaves the sequence unchanged."""
        seq = list("ACGTACGT")
        expected = list(seq)
        applyAmbiguities(seq, [])
        self.assertEqual(expected, seq)

    def testSingleHole(self):
        """A single N-run replaces the correct positions."""
        seq = list("AAAAAAAAAA")
        offset, length, char = 2, 3, "N"
        expected = list(seq)
        expected[offset : offset + length] = [char] * length
        applyAmbiguities(seq, [(offset, length, char)])
        self.assertEqual(expected, seq)

    def testHoleAtStart(self):
        """An N-run starting at offset 0 replaces positions from the beginning."""
        seq = list("ACGTACGT")
        offset, length, char = 0, 2, "N"
        expected = list(seq)
        expected[offset : offset + length] = [char] * length
        applyAmbiguities(seq, [(offset, length, char)])
        self.assertEqual(expected, seq)

    def testHoleAtEnd(self):
        """An N-run at the very end of the sequence replaces the final positions."""
        seq = list("ACGTACGT")
        offset, length, char = 6, 2, "N"
        expected = list(seq)
        expected[offset : offset + length] = [char] * length
        applyAmbiguities(seq, [(offset, length, char)])
        self.assertEqual(expected, seq)

    def testMultipleHoles(self):
        """Multiple N-runs all replace their respective positions."""
        seq = list("AAAAAAAAAA")
        holes = [(0, 1, "N"), (5, 2, "N")]
        expected = list(seq)
        for offset, length, char in holes:
            expected[offset : offset + length] = [char] * length
        applyAmbiguities(seq, holes)
        self.assertEqual(expected, seq)

    def testHoleDoesNotExceedSequenceLength(self):
        """A hole that extends beyond the sequence boundary does not raise an error."""
        seq = list("AAAA")
        offset, length, char = 2, 5, "N"
        # offset + length would reach index 7 which is out of range — the
        # guard `if pos < len(seq)` in applyAmbiguities handles this silently.
        expected = list(seq)
        expected[offset:] = [char] * (len(seq) - offset)
        applyAmbiguities(seq, [(offset, length, char)])
        self.assertEqual(expected, seq)

    def testNonNCharacter(self):
        """Ambiguity characters other than N are also applied correctly."""
        seq = list("ACGT")
        offset, length, char = 1, 2, "X"
        expected = list(seq)
        expected[offset : offset + length] = [char] * length
        applyAmbiguities(seq, [(offset, length, char)])
        self.assertEqual(expected, seq)


class TestReconstructFasta(TestCase):
    """Tests for reconstructFasta()."""

    def testSingleSequenceNoAnnotation(self):
        """A single sequence without an annotation writes the expected FASTA."""
        name = "chr1"
        bases = "ACGT"
        seqs = [{"name": name, "anno": "", "offset": 0, "length": len(bases)}]
        seqData = list(bases)
        result = reconstructFasta(seqs, [], seqData, lineWidth=60)
        expected = f">{name}\n{bases}\n"
        self.assertEqual(expected, result)

    def testSingleSequenceWithAnnotation(self):
        """An annotation is appended to the header line after a space."""
        name = "seq1"
        anno = "chromosome 1"
        bases = "ACGT"
        seqs = [{"name": name, "anno": anno, "offset": 0, "length": len(bases)}]
        seqData = list(bases)
        result = reconstructFasta(seqs, [], seqData)
        expected = f">{name} {anno}\n{bases}\n"
        self.assertEqual(expected, result)

    def testMultipleSequences(self):
        """Multiple sequences are all written in order."""
        name1, bases1 = "s1", "ACG"
        name2, bases2 = "s2", "TT"
        seqs = [
            {"name": name1, "anno": "", "offset": 0, "length": len(bases1)},
            {"name": name2, "anno": "", "offset": len(bases1), "length": len(bases2)},
        ]
        seqData = list(bases1 + bases2)
        result = reconstructFasta(seqs, [], seqData)
        expected = f">{name1}\n{bases1}\n>{name2}\n{bases2}\n"
        self.assertEqual(expected, result)

    def testLineWrapping(self):
        """Sequences longer than lineWidth are split across multiple lines."""
        name = "seq"
        bases = "ACGTACGTAC"
        lineWidth = 4
        seqs = [{"name": name, "anno": "", "offset": 0, "length": len(bases)}]
        seqData = list(bases)
        result = reconstructFasta(seqs, [], seqData, lineWidth=lineWidth)
        expected = f">{name}\nACGT\nACGT\nAC\n"
        self.assertEqual(expected, result)

    def testOffsetSlicesCorrectly(self):
        """The offset field in the sequence descriptor is honoured."""
        name = "s2"
        prefix = "AAAA"
        bases = "CCCC"
        seqs = [{"name": name, "anno": "", "offset": len(prefix), "length": len(bases)}]
        seqData = list(prefix + bases)
        result = reconstructFasta(seqs, [], seqData)
        expected = f">{name}\n{bases}\n"
        self.assertEqual(expected, result)

    def testEmptySequence(self):
        """A sequence of length 0 produces a header line but no bases."""
        name = "empty"
        seqs = [{"name": name, "anno": "", "offset": 0, "length": 0}]
        result = reconstructFasta(seqs, [], [])
        expected = f">{name}\n"
        self.assertEqual(expected, result)

    def testAmbiguousBasesInOutput(self):
        """N bases stored in seqData are written out as-is."""
        name = "s"
        bases = "ACNNGT"
        seqs = [{"name": name, "anno": "", "offset": 0, "length": len(bases)}]
        seqData = list(bases)
        result = reconstructFasta(seqs, [], seqData)
        expected = f">{name}\n{bases}\n"
        self.assertEqual(expected, result)


@skipUnless(BWA_AVAILABLE, "bwa is not available — skipping integration tests")
class TestBwaIntegration(TestCase):
    """Round-trip tests: build a real BWA index then reconstruct via the module."""

    def _runBwaIndex(self, faPath):
        """Run bwa index on *faPath*, raising on failure."""
        subprocess.run(
            ["bwa", "index", str(faPath)],
            check=True,
            capture_output=True,
        )

    def _reconstruct(self, faPath, lineWidth):
        """Read index files and reconstruct FASTA."""
        prefix = str(faPath)
        _, sequences = readAnn(f"{prefix}.ann")
        holes = readAmb(f"{prefix}.amb")
        seqData = readPac(f"{prefix}.pac")
        return reconstructFasta(sequences, holes, seqData, lineWidth)

    def testSimpleSequenceRoundTrip(self):
        """A plain sequence with no Ns is reconstructed exactly."""
        name = "chr1"
        bases = "ACGTACGTACGT"
        fasta = f">{name}\n{bases}\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = Path(tmpdir) / "ref.fa"
            fa.write_text(fasta)
            self._runBwaIndex(fa)
            reconstructed = self._reconstruct(fa, lineWidth=len(bases))

        self.assertEqual(fasta, reconstructed)

    def testSequenceWithNsRoundTrip(self):
        """Ambiguous N bases are restored to their original positions."""
        name = "chr1"
        bases = "ACGTNNACGT"
        fasta = f">{name}\n{bases}\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = Path(tmpdir) / "ref.fa"
            fa.write_text(fasta)
            self._runBwaIndex(fa)
            reconstructed = self._reconstruct(fa, lineWidth=len(bases))

        self.assertEqual(fasta, reconstructed)

    def testMultipleSequencesRoundTrip(self):
        """All sequences in a multi-record FASTA are reconstructed correctly."""
        sequences = [("seq1", "AAAA"), ("seq2", "CCCC"), ("seq3", "GGGG")]
        fasta = "".join(f">{name}\n{bases}\n" for name, bases in sequences)
        maxLen = max(len(bases) for _, bases in sequences)
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = Path(tmpdir) / "ref.fa"
            fa.write_text(fasta)
            self._runBwaIndex(fa)
            reconstructed = self._reconstruct(fa, lineWidth=maxLen)

        self.assertEqual(fasta, reconstructed)

    def testSequenceNamesPreserved(self):
        """Sequence names from the FASTA header are preserved in the output."""
        name = "my_special_contig"
        bases = "ACGTACGT"
        fasta = f">{name}\n{bases}\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = Path(tmpdir) / "ref.fa"
            fa.write_text(fasta)
            self._runBwaIndex(fa)
            reconstructed = self._reconstruct(fa, lineWidth=len(bases))

        self.assertEqual(fasta, reconstructed)

    def testLongerSequenceRoundTrip(self):
        """A sequence longer than one output line is reconstructed in full."""
        name = "long"
        bases = "ACGT" * 20  # 80 bases
        fasta = f">{name}\n{bases}\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = Path(tmpdir) / "ref.fa"
            fa.write_text(fasta)
            self._runBwaIndex(fa)
            reconstructed = self._reconstruct(fa, lineWidth=len(bases))

        self.assertEqual(fasta, reconstructed)

    def testMultipleNRunsRoundTrip(self):
        """Several separate N-runs in a single sequence are all restored."""
        name = "chr1"
        bases = "ACGTNNACGTNNNNACGT"
        fasta = f">{name}\n{bases}\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = Path(tmpdir) / "ref.fa"
            fa.write_text(fasta)
            self._runBwaIndex(fa)
            reconstructed = self._reconstruct(fa, lineWidth=len(bases))

        self.assertEqual(fasta, reconstructed)

    def testMixedSequencesWithNsRoundTrip(self):
        """Multiple sequences where only some contain Ns are all correct."""
        sequences = [("clean", "AAAACCCC"), ("ambig", "GGGNNNTTT")]
        fasta = "".join(f">{name}\n{bases}\n" for name, bases in sequences)
        maxLen = max(len(bases) for _, bases in sequences)
        with tempfile.TemporaryDirectory() as tmpdir:
            fa = Path(tmpdir) / "ref.fa"
            fa.write_text(fasta)
            self._runBwaIndex(fa)
            reconstructed = self._reconstruct(fa, lineWidth=maxLen)

        self.assertEqual(fasta, reconstructed)

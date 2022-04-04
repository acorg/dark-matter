from unittest import TestCase

from dark.aligners import edlibAlign
from dark.reads import DNARead, Reads


class TestEdlibAlign(TestCase):
    """
    Test the edlibAlign function.
    """
    def testEmtpyStrings(self):
        """
        Aligning two empty reads results in an Exception.
        """
        in1 = DNARead('id1', '')
        in2 = DNARead('id2', '')

        error = (r"^The object alignResult contains an empty CIGAR string\. "
                 r"Users must run align\(\) with task='path'\. Please check "
                 r"the input alignResult\.$")
        self.assertRaisesRegex(Exception, error, edlibAlign, Reads([in1, in2]))

    def testIdenticalStringsOfLengthOne(self):
        """
        Aligning identical reads of length one must produce the expected
        result.
        """
        in1 = DNARead('id1', 'A')
        in2 = DNARead('id2', 'A')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testNonIdenticalStringsOfLengthOne(self):
        """
        Aligning non-identical reads of length one must produce the expected
        result.
        """
        in1 = DNARead('id1', 'A')
        in2 = DNARead('id2', 'G')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testIdenticalStringsOfLengthTwo(self):
        """
        Aligning identical reads of length two must produce the expected
        result.
        """
        in1 = DNARead('id1', 'AG')
        in2 = DNARead('id2', 'AG')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testNonIdenticalStringsOfLengthTwo(self):
        """
        Aligning non-identical reads of length two must produce the expected
        result.
        """
        in1 = DNARead('id1', 'AT')
        in2 = DNARead('id2', 'GC')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

    def testGappedPrefix(self):
        """
        Aligning a read that is a suffix of another must result in gaps being
        placed at the start of the shorter string.
        """
        in1 = DNARead('id1', 'ATGCCGTTCA')
        in2 = DNARead('id2', '--GCCGTTCA')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual('--GCCGTTCA', out2.sequence)

    def testLongerGappedPrefix(self):
        """
        Aligning a read that is a suffix of another must result in gaps being
        placed at the expected places at the start of the shorter string.
        """
        in1 = DNARead('id1', 'ATGCCGTTCA')
        in2 = DNARead('id2', 'CGTTCA')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual('---C-GTTCA', out2.sequence)

    def testGappedSuffix(self):
        """
        Aligning a read that is a prefix of another must result in gaps being
        placed at the end of the shorter string.
        """
        in1 = DNARead('id1', 'ATGCCGTTCA')
        in2 = DNARead('id2', 'ATGCCGTT')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual('ATGCCGTT--', out2.sequence)

    def testLongerGappedSuffix(self):
        """
        Aligning a read that is a prefix of another must result in gaps being
        placed at the expected places at the end of the shorter string.
        """
        in1 = DNARead('id1', 'ATGCCGTTCA')
        in2 = DNARead('id2', 'ATGCC')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual('ATGCC-----', out2.sequence)

    def testLongerGappedSuffixWithGap(self):
        """
        Aligning a read that is a prefix of another must result in gaps being
        placed at the expected places at the end of the shorter string.
        """
        in1 = DNARead('id1', 'ATGCCGTTCA')
        in2 = DNARead('id2', 'ATGCCTT')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual('ATGCC-TT--', out2.sequence)

    def testGapsAtBothEnds(self):
        """
        Aligning two reads that have a common prefix/suffix must result
        in gaps at the start of one and the end of the other.
        """
        in1 = DNARead('id1', 'XXXATGCCGTTCA')
        in2 = DNARead('id2', 'ATGCCGTTCAYYY')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual('XXXATGCCGTTCA---', out1.sequence)
        self.assertEqual('---ATGCCGTTCAYYY', out2.sequence)

    def testGapsAtBothEndsAndOneInTheMiddle(self):
        """
        Aligning two reads that have a common prefix/suffix and one deletion
        must result in gaps at the start of one and the end of the other,
        and one in the middle.
        """
        in1 = DNARead('id1', 'XXXATGCCTTCA')
        in2 = DNARead('id2', 'ATGCCGTTCAYYY')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual('XXXATGCC-TTCA---', out1.sequence)
        self.assertEqual('---ATGCCGTTCAYYY', out2.sequence)

    def testGapsAtBothEndsAndTwoInTheMiddle(self):
        """
        Aligning two reads that have a common prefix/suffix and two deletions
        must result in gaps at the start of one and the end of the other,
        and two in the middle.
        """
        in1 = DNARead('id1', 'XXXATCCTTCA')
        in2 = DNARead('id2', 'ATGCCGTTCAYYY')

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual('XXXAT-CC-TTCA---', out1.sequence)
        self.assertEqual('---ATGCCGTTCAYYY', out2.sequence)

    def test30Knucleotides(self):
        """
        Aligning two reads of length 30K nucleotides must work (and run
        quickly).
        """
        in1 = DNARead('id1', 'A' * 30_000)
        in2 = DNARead('id2', 'A' * 30_000)

        out1, out2 = list(edlibAlign(Reads([in1, in2])))

        self.assertEqual(in1, out1)
        self.assertEqual(in2, out2)

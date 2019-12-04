import six
from unittest import TestCase

from dark.reads import Read, Reads
from dark.score import HigherIsBetterScore
from dark.hsp import HSP, LSP
from dark.alignments import (
    Alignment, bestAlignment, ReadAlignments, ReadsAlignmentsParams,
    ReadsAlignments)


class TestAlignment(TestCase):
    """
    Tests for the dark.alignment.Alignment class
    """

    def testExpectedAttrs(self):
        """
        An alignment must have the expected attributes.
        """
        alignment = Alignment(45, 'title')
        self.assertEqual('title', alignment.subjectTitle)
        self.assertEqual(45, alignment.subjectLength)

    def testNoHspsWhenCreated(self):
        """
        An alignment must have no HSPs when it is created.
        """
        alignment = Alignment(45, 'title')
        self.assertEqual(0, len(alignment.hsps))

    def testAddHsp(self):
        """
        It must be possible to add an HSP to an alignment.
        """
        alignment = Alignment(45, 'title')
        alignment.addHsp(HSP(3))
        self.assertEqual(HSP(3), alignment.hsps[0])


class TestReadAlignments(TestCase):
    """
    Tests for the dark.alignment.ReadAlignments class
    """

    def testRead(self):
        """
        An read alignments must store its read.
        """
        read = Read('id', 'ACGT')
        readAlignments = ReadAlignments(read)
        self.assertEqual(read, readAlignments.read)

    def testNoAlignments(self):
        """
        An read alignments must be able to have no alignments.
        """
        read = Read('id', 'ACGT')
        readAlignments = ReadAlignments(read)
        self.assertEqual(0, len(readAlignments))

    def testAlignments(self):
        """
        An read alignments must store its alignments.
        """
        read = Read('id', 'ACGT')
        alignment1 = Alignment(45, 'title1')
        alignment2 = Alignment(55, 'title2')
        readAlignments = ReadAlignments(read, [alignment1, alignment2])
        self.assertEqual([alignment1, alignment2], readAlignments)


class TestBestAlignmentHSP(TestCase):
    """
    Test the L{dark.hits.bestAlignment} function when HSPs are used.
    """

    def testOneAlignment(self):
        """
        When one alignment is present that alignment must be returned by
        bestAlignment.
        """
        alignment = Alignment(44, 'Seq 1')
        alignment.addHsp(HSP(10))
        alignment.addHsp(HSP(9))

        alignments = [alignment]
        hit = ReadAlignments(Read('id1', 'aaa'), alignments)
        best = bestAlignment(hit)
        self.assertEqual('Seq 1', best.subjectTitle)
        self.assertEqual(44, best.subjectLength)

    def testThreeAlignments(self):
        """
        When three alignments are present, the one with the highest first HSP
        must be returned by bestAlignment.
        """
        alignment1 = Alignment(33, 'Seq 1')
        alignment1.addHsp(HSP(10))
        alignment1.addHsp(HSP(9))

        alignment2 = Alignment(44, 'Seq 2')
        alignment2.addHsp(HSP(30))
        alignment2.addHsp(HSP(29))

        alignment3 = Alignment(55, 'Seq 3')
        alignment3.addHsp(HSP(20))
        alignment3.addHsp(HSP(19))

        alignments = [alignment1, alignment2, alignment3]
        hit = ReadAlignments(Read('id1', 'aaa'), alignments)
        best = bestAlignment(hit)
        self.assertEqual('Seq 2', best.subjectTitle)
        self.assertEqual(44, best.subjectLength)


class TestBestAlignmentLSP(TestCase):
    """
    Test the L{dark.hits.bestAlignment} function when LSPs are used.
    """

    def testOneAlignment(self):
        """
        When one alignment is present that alignment must be returned by
        bestAlignment.
        """
        alignment = Alignment(44, 'Seq 1')
        alignment.addHsp(LSP(10))
        alignment.addHsp(LSP(9))

        alignments = [alignment]
        readAlignments = ReadAlignments(Read('id0', 'aaa'), alignments)
        best = bestAlignment(readAlignments)
        self.assertEqual('Seq 1', best.subjectTitle)
        self.assertEqual(44, best.subjectLength)

    def testThreeAlignments(self):
        """
        When three alignments are present, the one with the lowest first HSP
        must be returned by bestAlignment.
        """
        alignment1 = Alignment(33, 'Seq 1')
        alignment1.addHsp(LSP(10))
        alignment1.addHsp(LSP(9))

        alignment2 = Alignment(44, 'Seq 2')
        alignment2.addHsp(LSP(3))
        alignment2.addHsp(LSP(2))

        alignment3 = Alignment(55, 'Seq 3')
        alignment3.addHsp(LSP(20))
        alignment3.addHsp(LSP(19))

        alignments = [alignment1, alignment2, alignment3]
        readAlignments = ReadAlignments(Read('id0', 'aaa'), alignments)
        best = bestAlignment(readAlignments)
        self.assertEqual('Seq 2', best.subjectTitle)
        self.assertEqual(44, best.subjectLength)


class TestReadsAlignmentsParams(TestCase):
    """
    Test the L{dark.alignments.ReadsAlignmentsParams} class.
    """

    def testExpectedAttrs(self):
        """
        A ReadsAlignmentsParams instance must have the expected attributes.
        """
        applicationParams = {}
        params = ReadsAlignmentsParams('application name', applicationParams,
                                       False, 'Bit score')
        self.assertEqual('application name', params.application)
        self.assertIs(applicationParams, params.applicationParams)
        self.assertFalse(params.subjectIsNucleotides)
        self.assertEqual('Bit score', params.scoreTitle)


class TestReadsAlignments(TestCase):
    """
    Test the L{dark.alignments.ReadsAlignments} class.
    """

    # NOTE: The ReadsAlignments class is a base class for concrete
    # implementations, such as BlastReadsAlignments. So it can only be
    # tested minimally by itself. For full(er) tests see
    # diamond/test_alignments.py and blast/test_alignments.py.

    def testExpectedAttrs(self):
        """
        A ReadsAlignments instance must have the expected attributes.
        """
        reads = Reads()
        params = {
            'application': 'app name'
        }
        readsAlignments = ReadsAlignments(reads, params)
        self.assertIs(readsAlignments.reads, reads)
        self.assertEqual('app name', readsAlignments.params['application'])
        self.assertIs(params, readsAlignments.params)
        self.assertIs(HigherIsBetterScore, readsAlignments.scoreClass)

    def testNotIterable(self):
        """
        Iterating an empty ReadsAlignments must result in the empty list.
        """
        reads = Reads()
        readsAlignments = ReadsAlignments(reads, 'applicationName', None)
        self.assertEqual([], list(readsAlignments))

    def testGetSubjectSequence(self):
        """
        A ReadsAlignments instance will not implement getSubjectSequence.
        Subclasses are expected to implement it.
        """
        reads = Reads()
        readsAlignments = ReadsAlignments(reads, 'applicationName', None)
        error = 'getSubjectSequence must be implemented by a subclass'
        six.assertRaisesRegex(self, NotImplementedError, error,
                              readsAlignments.getSubjectSequence, 'title')

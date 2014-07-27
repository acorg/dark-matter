from unittest import TestCase

from dark.reads import Read, Reads
from dark.hsp import HSP
from dark.alignments import (
    Alignment, bestAlignment, ReadAlignments, ReadsAlignments)


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
        # TODO: make this a proper HSP.
        alignment.addHsp(3)
        self.assertEqual(3, alignment.hsps[0])


class TestBestAlignment(TestCase):
    """
    Test the L{dark.hits.bestAlignment} function.
    """

    def testOneAlignment(self):
        """
        When one alignment is present that alignment must be returned by
        bestAlignment.
        """
        alignment = Alignment(44, 'Seq 1')
        alignment.addHsp(HSP(score=10))
        alignment.addHsp(HSP(score=9))

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
        alignment1.addHsp(HSP(score=10))
        alignment1.addHsp(HSP(score=9))

        alignment2 = Alignment(44, 'Seq 2')
        alignment2.addHsp(HSP(score=30))
        alignment2.addHsp(HSP(score=29))

        alignment3 = Alignment(55, 'Seq 3')
        alignment3.addHsp(HSP(score=20))
        alignment3.addHsp(HSP(score=19))

        alignments = [alignment1, alignment2, alignment3]
        hit = ReadAlignments(Read('id1', 'aaa'), alignments)
        best = bestAlignment(hit)
        self.assertEqual('Seq 2', best.subjectTitle)
        self.assertEqual(44, best.subjectLength)


class TestReadsAlignments(TestCase):
    """
    Test the L{dark.hits.ReadsAlignments} class.
    """

    def testNotIterable(self):
        """
        A ReadsAlignments instance must not be filterable. A subclass is
        expected to implement __iter__.
        """
        reads = Reads()
        hits = ReadsAlignments(reads, 'applicationName', None)
        error = '__iter__ must be implemented by a subclass'
        self.assertRaisesRegexp(NotImplementedError, error, list,
                                hits.filter())

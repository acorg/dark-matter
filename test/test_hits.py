from unittest import TestCase

from dark.reads import Read
from dark.hsp import HSP
from dark.hits import Hit, Hits, bestAlignment
from dark.alignment import Alignment


class TestBestAlignment(TestCase):
    """
    Test the bestAlignment function.
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
        hit = Hit(Read('id1', 'aaa'), alignments)
        best = bestAlignment(hit)
        self.assertEqual('Seq 1', best.hitTitle)
        self.assertEqual(44, best.hitLength)

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
        hit = Hit(Read('id1', 'aaa'), alignments)
        best = bestAlignment(hit)
        self.assertEqual('Seq 2', best.hitTitle)
        self.assertEqual(44, best.hitLength)


class TestReadHits(TestCase):
    """
    Test the ReadHits class.
    """

    def testNoHits(self):
        """
        A ReadHits instance must not be filterable. A subclass is expected
        to implement __iter__.
        """
        # readHits = ReadHits()
        # self.assertRaises(TypeError, list, readHits.filter())

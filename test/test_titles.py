from unittest import TestCase

from dark.titles import TitleAlignment, TitleAlignments
from dark.reads import Read
from dark.hsp import HSP


class TestTitleAlignment(TestCase):
    """
    Test the TitleAlignment class
    """

    def testExpectedAttributes(self):
        """
        An instance of TitleAlignment must have the expected attributes.
        """
        read = Read('id', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        self.assertEqual(read, titleAlignment.read)
        self.assertEqual([], titleAlignment.hsps)


class TestTitleAlignments(TestCase):
    """
    Test the TitleAlignments class
    """

    def testExpectedAttributes(self):
        """
        An instance of TitleAlignments must have the expected attributes.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        self.assertEqual('subject title', titleAlignments.subjectTitle)
        self.assertEqual(55, titleAlignments.subjectLength)
        self.assertEqual([], titleAlignments.alignments)

    def testAddAlignment(self):
        """
        It must be possible to add an alignment to an instance of
        TitleAlignments.
        """
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id', 'AAA')
        titleAlignment = TitleAlignment(read, [])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual(read, titleAlignments.alignments[0].read)
        self.assertEqual([], titleAlignments.alignments[0].hsps)

    def testHSPs(self):
        """
        The hsps function must produce a list of all HSPs.
        """
        hsp1 = HSP(readStart=1, readEnd=2,
                   readStartInHit=3, readEndInHit=4,
                   hitStart=5, hitEnd=6,
                   readMatchedSequence='aaa', hitMatchedSequence='ccc',
                   score=7)
        hsp2 = HSP(readStart=8, readEnd=9,
                   readStartInHit=10, readEndInHit=11,
                   hitStart=12, hitEnd=13,
                   readMatchedSequence='aaa', hitMatchedSequence='ccc',
                   score=14)
        hsp3 = HSP(readStart=15, readEnd=16,
                   readStartInHit=17, readEndInHit=18,
                   hitStart=19, hitEnd=20,
                   readMatchedSequence='aaa', hitMatchedSequence='ccc',
                   score=21)
        titleAlignments = TitleAlignments('subject title', 55)
        read = Read('id1', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp1, hsp2])
        titleAlignments.addAlignment(titleAlignment)
        read = Read('id2', 'AAA')
        titleAlignment = TitleAlignment(read, [hsp3])
        titleAlignments.addAlignment(titleAlignment)
        self.assertEqual([7, 14, 21],
                         [hsp.score for hsp in titleAlignments.hsps()])

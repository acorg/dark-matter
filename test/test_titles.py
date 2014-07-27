from unittest import TestCase

from dark.titles import TitleAlignment, TitleAlignments
from dark.reads import Read


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

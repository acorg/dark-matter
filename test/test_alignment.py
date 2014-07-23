from unittest import TestCase

from dark.alignment import Alignment


class TestAlignment(TestCase):
    """
    Tests for the dark.alignment.Alignment class
    """

    def testExpectedAttrs(self):
        """
        An alignment must have the expected attributes.
        """
        alignment = Alignment(45, 'title')
        self.assertEqual('title', alignment.hitTitle)
        self.assertEqual(45, alignment.hitLength)

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

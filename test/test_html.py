from unittest import TestCase
from six import assertRaisesRegex

from dark.html import NCBISequenceLinkURL, NCBISequenceLink


class TestNCBISequenceLinkURL(TestCase):
    """
    Test the NCBISequenceLinkURL function.
    """

    def testNoField(self):
        """
        If no field is passed, the passed title must be used.
        """
        self.assertEqual('http://www.ncbi.nlm.nih.gov/nuccore/xxx',
                         NCBISequenceLinkURL('xxx'))

    def testFieldNumber(self):
        title = 'gi|323924|gb|M15204.1|FCVMYCCA Feline leukemia virus myc gene'
        self.assertEqual('http://www.ncbi.nlm.nih.gov/nuccore/M15204.1',
                         NCBISequenceLinkURL(title, 3))

    def testAlternateDelimiter(self):
        title = 'gi+37955203+gb+AY253278.1+ Homo sapiens clone AL-11 HIV-1'
        self.assertEqual('http://www.ncbi.nlm.nih.gov/nuccore/AY253278.1',
                         NCBISequenceLinkURL(title, 3, '+'))

    def testFieldOutOfRange(self):
        """
        If the field number passed does not exist, IndexError must be raised.
        """
        title = 'gi|37955203|gb|AY253278.1| Homo sapiens clone AL-11 HIV-1'
        error = (
            r"^Could not extract field 10 from sequence title "
            r"'gi\|37955203\|gb\|AY253278\.1\| Homo sapiens clone "
            "AL-11 HIV-1'$")
        assertRaisesRegex(self, IndexError, error, NCBISequenceLinkURL, title,
                          10)


class TestNCBISequenceLink(TestCase):
    """
    Test the NCBISequenceLink function.
    """
    def testNoField(self):
        """
        If no field is passed, the passed title must be used.
        """
        self.assertEqual(
            '<a href="http://www.ncbi.nlm.nih.gov/nuccore/xxx" ' +
            'target="_blank">xxx</a>', NCBISequenceLink('xxx'))

    def testFieldNumber(self):
        title = 'gi|323924|gb|M15204.1|FCVMYCCA Feline leukemia virus myc gene'
        self.assertEqual(
            '<a href="http://www.ncbi.nlm.nih.gov/nuccore/M15204.1" ' +
            'target="_blank">' + title + '</a>', NCBISequenceLink(title, 3))

    def testAlternateDelimiter(self):
        title = 'gi+37955203+gb+AY253278.1+ Homo sapiens clone AL-11 HIV-1'
        self.assertEqual(
            '<a href="http://www.ncbi.nlm.nih.gov/nuccore/AY253278.1" ' +
            'target="_blank">' + title + '</a>',
            NCBISequenceLink(title, 3, '+'))

    def testFieldOutOfRange(self):
        """
        If the field number passed does not exist, IndexError must be raised.
        """
        title = 'gi|37955203|gb|AY253278.1| Homo sapiens clone AL-11 HIV-1'
        error = (
            r"^Could not extract field 10 from sequence title "
            r"'gi\|37955203\|gb\|AY253278\.1\| Homo sapiens clone "
            "AL-11 HIV-1'$")
        assertRaisesRegex(self, IndexError, error, NCBISequenceLink, title, 10)

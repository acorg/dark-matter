from unittest import TestCase
from dark.html import NCBISequenceLinkURL, NCBISequenceLink


class TestNCBISequenceLinkURL(TestCase):
    """
    Test the NCBISequenceLinkURL function.
    """

    def testGenericTitle1(self):
        title = 'gi|323924|gb|M15204.1|FCVMYCCA Feline leukemia virus myc gene'
        self.assertEqual('http://www.ncbi.nlm.nih.gov/nuccore/M15204',
                         NCBISequenceLinkURL(title))

    def testGenericTitle2(self):
        title = 'gi|37955203|gb|AY253278.1| Homo sapiens clone AL-11 HIV-1'
        self.assertEqual('http://www.ncbi.nlm.nih.gov/nuccore/AY253278',
                         NCBISequenceLinkURL(title))


class TestNCBISequenceLink(TestCase):
    """
    Test the NCBISequenceLink function.
    """

    def testGenericTitle1(self):
        title = 'gi|323924|gb|M15204.1|FCVMYCCA Feline leukemia virus myc gene'
        self.assertEqual(
            '<a href="http://www.ncbi.nlm.nih.gov/nuccore/M15204" ' +
            'target="_blank">' + title + '</a>', NCBISequenceLink(title))

    def testGenericTitle2(self):
        title = 'gi|37955203|gb|AY253278.1| Homo sapiens clone AL-11 HIV-1'
        self.assertEqual(
            '<a href="http://www.ncbi.nlm.nih.gov/nuccore/AY253278" ' +
            'target="_blank">' + title + '</a>', NCBISequenceLink(title))

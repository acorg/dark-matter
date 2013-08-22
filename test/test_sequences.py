from unittest import TestCase
from mock import patch

from mocking import mockOpen
from dark.sequences import summarizeRecordsBySequence


class TestSummarizeRecordsBySequence(TestCase):
    """
    Test the summarizeRecordsBySequence function.
    """

    def testEmptyJSONInput(self):
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            self.assertEqual({}, summarizeRecordsBySequence('f.json'))

    def testEmptyXMLInput(self):
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            self.assertRaises(ValueError, summarizeRecordsBySequence, 'f.xml')

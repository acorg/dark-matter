from json import dumps
from unittest import TestCase
from mock import patch

from ..mocking import mockOpen
from sample_data import PARAMS, RECORD0, RECORD1, RECORD2, RECORD3

from dark.reads import Read, Reads
from dark.blast.alignments import BlastReadsAlignments
from dark.titles import titleCounts, TitlesAlignments


class TestTitleCounts(TestCase):
    """
    Test the titleCounts function.
    """

    def testEmpty(self):
        """
        If passed an empty readsAlignments, titleCounts must return an
        empty dictionary.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual({}, titleCounts(readsAlignments))

    def testThreeRecords(self):
        """
        If alignments for three reads are passed to titleCounts, it must
        return the correct title counts.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual(
                {
                    'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99': 1,
                    'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.': 1,
                    'gi|887699|gb|DQ37780 Cowpox virus 15': 1,
                    'gi|887699|gb|DQ37780 Monkeypox virus 456': 1,
                    'gi|887699|gb|DQ37780 Squirrelpox virus 55': 1
                },
                titleCounts(readsAlignments))

    def testDuplicatedTitle(self):
        """
        If alignments for reads have a common title, the count on that title
        must be correct.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual(
                {
                    'gi|887699|gb|DQ37780 Cowpox virus 15': 2,
                },
                titleCounts(readsAlignments))


class TestTitlesAlignments(TestCase):
    """
    Test the TitlesAlignments class
    """

    def testEmpty(self):
        """
        An instance of TitlesAlignments must have no titles if passed an
        empty readsAlignments instance.
        """
        mockOpener = mockOpen(read_data=(dumps(PARAMS) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            self.assertEqual([], titlesAlignments.keys())

    def testExpectedTitles(self):
        """
        An instance of TitlesAlignments must have the expected titles.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Cowpox virus 15',
                    'gi|887699|gb|DQ37780 Monkeypox virus 456',
                    'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',
                    'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                    'gi|887699|gb|DQ37780 Squirrelpox virus 55',
                ],
                sorted(titlesAlignments.keys()))

    def testExpectedTitleDetails(self):
        """
        An instance of TitleAlignments in a TitlesAlignments instance must
        have the expected attributes.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            read = Read('id0', 'A' * 70)
            reads.add(read)
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)

            title = 'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99'
            titleAlignments = titlesAlignments[title]
            self.assertEqual(title, titleAlignments.subjectTitle)
            self.assertEqual(37000, titleAlignments.subjectLength)
            self.assertEqual(1, len(titleAlignments.alignments))
            self.assertEqual(read, titleAlignments.alignments[0].read)
            self.assertEqual(20, titleAlignments.alignments[0].hsps[0].score)

            title = 'gi|887699|gb|DQ37780 Squirrelpox virus 55'
            titleAlignments = titlesAlignments[title]
            self.assertEqual(title, titleAlignments.subjectTitle)
            self.assertEqual(38000, titleAlignments.subjectLength)
            self.assertEqual(1, len(titleAlignments.alignments))
            self.assertEqual(read, titleAlignments.alignments[0].read)
            self.assertEqual(25, titleAlignments.alignments[0].hsps[0].score)

    def testTitleCollection(self):
        """
        A title that occurs in the alignments of multiple reads must have
        the data from both reads collected properly.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            read2 = Read('id2', 'A' * 70)
            read3 = Read('id3', 'A' * 70)
            reads.add(read2)
            reads.add(read3)
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)

            title = 'gi|887699|gb|DQ37780 Cowpox virus 15'
            titleAlignments = titlesAlignments[title]
            self.assertEqual(title, titleAlignments.subjectTitle)
            self.assertEqual(30000, titleAlignments.subjectLength)
            self.assertEqual(2, len(titleAlignments.alignments))

            self.assertEqual(read2, titleAlignments.alignments[0].read)
            self.assertEqual(20, titleAlignments.alignments[0].hsps[0].score)

            self.assertEqual(read3, titleAlignments.alignments[1].read)
            self.assertEqual(20, titleAlignments.alignments[1].hsps[0].score)

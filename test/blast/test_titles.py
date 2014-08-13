from json import dumps
from unittest import TestCase
from mock import patch

from ..mocking import mockOpen
from sample_data import PARAMS, RECORD0, RECORD1, RECORD2, RECORD3, RECORD4

from dark.reads import Read, Reads
from dark.hsp import HSP
from dark.score import LowerIsBetterScore
from dark.blast.alignments import BlastReadsAlignments
from dark.titles import titleCounts, TitleAlignments, TitlesAlignments


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
            self.assertEqual(HSP(20), titleAlignments.alignments[0].hsps[0])

            title = 'gi|887699|gb|DQ37780 Squirrelpox virus 55'
            titleAlignments = titlesAlignments[title]
            self.assertEqual(title, titleAlignments.subjectTitle)
            self.assertEqual(38000, titleAlignments.subjectLength)
            self.assertEqual(1, len(titleAlignments.alignments))
            self.assertEqual(read, titleAlignments.alignments[0].read)
            self.assertEqual(HSP(25), titleAlignments.alignments[0].hsps[0])

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
            self.assertEqual(HSP(20), titleAlignments.alignments[0].hsps[0])

            self.assertEqual(read3, titleAlignments.alignments[1].read)
            self.assertEqual(HSP(20), titleAlignments.alignments[1].hsps[0])

    def testAddTitleRepeat(self):
        """
        The addTitle function must raise a C{KeyError} if an attempt is made
        to add a pre-existing title to a TitlesAlignments instance.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            title = 'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99'
            titleAlignments = TitleAlignments(title, 55)
            error = ("Title 'gi\|887699\|gb\|DQ37780 Squirrelpox virus "
                     "1296/99' already present in TitlesAlignments instance\.")
            self.assertRaisesRegexp(KeyError, error, titlesAlignments.addTitle,
                                    title, titleAlignments)

    def testAddTitle(self):
        """
        The addTitle function must add a title to the TitlesAlignments
        instance.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            title = 'gi|887699|gb|DQ37780 Squirrelpox virus 23'
            titleAlignments = TitleAlignments(title, 55)
            self.assertTrue(title not in titlesAlignments)
            titlesAlignments.addTitle(title, titleAlignments)
            self.assertTrue(title in titlesAlignments)

    def testHsps(self):
        """
        The hsps function must yield all the hsps for all titles in a
        TitlesAlignments instance.
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
            result = list(titlesAlignments.hsps())
            self.assertEqual(
                sorted([HSP(20), HSP(25), HSP(20), HSP(20), HSP(20)]),
                sorted(result))


class TestTitlesAlignmentsFiltering(TestCase):
    """
    Test the TitlesAlignments class filter function.
    """

    def testFilterWithNoArguments(self):
        """
        The filter function must return a TitlesAlignments instance with all
        the titles of the original when called with no arguments.
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
            result = titlesAlignments.filter()
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Cowpox virus 15',
                    'gi|887699|gb|DQ37780 Monkeypox virus 456',
                    'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',
                    'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                    'gi|887699|gb|DQ37780 Squirrelpox virus 55',
                ],
                sorted(result.keys()))

    def testMinMatchingReads(self):
        """
        The filter function work correctly when passed a value for
        minMatchingReads.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minMatchingReads=2)
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Cowpox virus 15',
                ],
                result.keys())

    def testMinMedianScore_Bits(self):
        """
        The filter function work correctly when passed a value for
        minMedianScore when using bit scores.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minMedianScore=22)
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Squirrelpox virus 55',
                ],
                result.keys())

    def testMinMedianScore_EValue(self):
        """
        The filter function work correctly when passed a value for
        minMedianScore when using e values.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minMedianScore=1e-9)
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                    'gi|887699|gb|DQ37780 Squirrelpox virus 55',
                ],
                sorted(result.keys()))

    def testWithScoreBetterThan_Bits(self):
        """
        The filter function work correctly when passed a value for
        withScoreBetterThan when using bit scores.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(withScoreBetterThan=24)
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Squirrelpox virus 55',
                ],
                result.keys())

    def testWithScoreBetterThan_EValue(self):
        """
        The filter function work correctly when passed a value for
        withScoreBetterThan when using e values.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(withScoreBetterThan=1e-10)
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                ],
                result.keys())

    def testReadSetFilterAllowAnything(self):
        """
        The filter function work correctly when passed a 0.0 value for
        minNewReads, i.e. that considers any read set sufficiently novel.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minNewReads=0.0)
            self.assertEqual(
                [
                    'gi|887699|gb|DQ37780 Cowpox virus 15',
                    'gi|887699|gb|DQ37780 Monkeypox virus 456',
                    'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',
                    'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                    'gi|887699|gb|DQ37780 Squirrelpox virus 55',
                ],
                sorted(result.keys()))

    def testReadSetFilterStrict(self):
        """
        The filter function work correctly when passed a 1.0 value for
        minNewReads.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.filter(minNewReads=1.0)

            # Either 'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.'
            # invalidates 'gi|887699|gb|DQ37780 Monkeypox virus 456' or
            # vice-versa. It depends on Python's dict walking order. Check
            # for both, making sure just one of them is true.

            mummypox = 'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.'
            monkeypox = 'gi|887699|gb|DQ37780 Monkeypox virus 456'

            assertionCount = 0
            if mummypox in result:
                self.assertTrue(monkeypox in
                                result.readSetFilter.invalidates(mummypox))
                assertionCount += 1
            if monkeypox in result:
                self.assertTrue(mummypox in
                                result.readSetFilter.invalidates(monkeypox))
                assertionCount += 1

            self.assertEqual(1, assertionCount)


class TestTitleSorting(TestCase):
    """
    Tests for the L{dark.titles.TitlesAlignments.sortTitles} function.
    """

    def testUnknown(self):
        """
        Sorting on an unknown attribute must raise C{ValueError}.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            self.assertRaises(ValueError, titlesAlignments.sortTitles, 'xxx')

    def testEmpty(self):
        """
        Sorting when there are no titles must return the empty list.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('title')
            self.assertEqual([], result)

    def testMedianScore_Bits(self):
        """
        Sorting on median score must work when scores are bit scores,
        including a secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n' + dumps(RECORD4) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            reads.add(Read('id4', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('medianScore')
            self.assertEqual([
                'gi|887699|gb|DQ37780 Squirrelpox virus 55',       # 25
                'gi|887699|gb|DQ37780 Monkeypox virus 456',        # 20
                'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',   # 20
                'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',  # 20
                'gi|887699|gb|DQ37780 Cowpox virus 15',            # 15
            ], result)

    def testMedianScore_EValue(self):
        """
        Sorting on median score must work when scores are bit scores,
        including a secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n' + dumps(RECORD4) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            reads.add(Read('id4', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('medianScore')
            self.assertEqual([
                'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',  # 1e-11
                'gi|887699|gb|DQ37780 Squirrelpox virus 55',       # 1e-10
                'gi|887699|gb|DQ37780 Monkeypox virus 456',        # 1e-8
                'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',   # 1e-7
                'gi|887699|gb|DQ37780 Cowpox virus 15',            # worst :-)
            ], result)

    def testMaxScore_Bits(self):
        """
        Sorting on max score must work when scores are bit scores, including a
        secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('maxScore')
            self.assertEqual([
                'gi|887699|gb|DQ37780 Squirrelpox virus 55',       # 25
                'gi|887699|gb|DQ37780 Cowpox virus 15',            # 20
                'gi|887699|gb|DQ37780 Monkeypox virus 456',        # 20
                'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',   # 20
                'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',  # 20
            ], result)

    def testMaxScore_EValue(self):
        """
        Sorting on max score must work when scores are e values, including a
        secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('maxScore')
            self.assertEqual([
                'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',  # 1e-11
                'gi|887699|gb|DQ37780 Squirrelpox virus 55',       # 1e-10
                'gi|887699|gb|DQ37780 Monkeypox virus 456',        # 1e-8
                'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',   # 1e-7
                'gi|887699|gb|DQ37780 Cowpox virus 15',            # 1e-6
            ], result)

    def testReadCount(self):
        """
        Sorting on read count must work, including a secondary sort on title.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('readCount')
            self.assertEqual([
                'gi|887699|gb|DQ37780 Cowpox virus 15',            # 3
                'gi|887699|gb|DQ37780 Monkeypox virus 456',        # 1
                'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',   # 1
                'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',  # 1
                'gi|887699|gb|DQ37780 Squirrelpox virus 55',       # 1
            ], result)

    def testLength(self):
        """
        Sorting on sequence length must work, including a secondary sort on
        title.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('length')
            self.assertEqual([
                'gi|887699|gb|DQ37780 Squirrelpox virus 55',       # 38000
                'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',  # 37000
                'gi|887699|gb|DQ37780 Monkeypox virus 456',        # 35000
                'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',   # 35000
                'gi|887699|gb|DQ37780 Cowpox virus 15',            # 30000
            ], result)

    def testTitle(self):
        """
        Sorting on title must work.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch('__builtin__.open', mockOpener, create=True):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            result = titlesAlignments.sortTitles('title')
            self.assertEqual([
                'gi|887699|gb|DQ37780 Cowpox virus 15',
                'gi|887699|gb|DQ37780 Monkeypox virus 456',
                'gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',
                'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                'gi|887699|gb|DQ37780 Squirrelpox virus 55',
            ], result)

# TODO: Add tests based on taxonomy, once we know how to mock mysql.

import six
import platform
from six.moves import builtins
from copy import deepcopy
from json import dumps
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from ..mocking import mockOpen, File
from .sample_data import PARAMS, RECORD0, RECORD1, RECORD2, RECORD3, RECORD4

from dark.reads import Read, Reads, AAReadWithX
from dark.hsp import HSP, LSP
from dark.score import LowerIsBetterScore
from dark.diamond.alignments import (
    DiamondReadsAlignments, ZERO_EVALUE_UPPER_RANDOM_INCREMENT)
from dark.titles import TitlesAlignments


class TestDiamondReadsAlignments(TestCase):
    """
    Test the DiamondReadsAlignments class.
    """

    def testEmptyJSONInput(self):
        """
        When a JSON input file is empty, a C{ValueError} must be raised
        on trying to read it.
        """
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            error = "JSON file 'file.json' was empty."
            six.assertRaisesRegex(self, ValueError, error,
                                  DiamondReadsAlignments, reads, 'file.json',
                                  'database.fasta')

    def testNonJSONInput(self):
        """
        When given a file whose contents are not JSON, attempting to
        read the DIAMOND hits from it must raise a C{ValueError}.
        """
        pypy = platform.python_implementation() == 'PyPy'
        mockOpener = mockOpen(read_data='not JSON\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            if six.PY3:
                error = (
                    "^Could not convert first line of 'file\.json' to JSON "
                    "\(Expecting value: line 1 column 1 \(char 0\)\)\. "
                    "Line is 'not JSON'\.$")
            else:
                if pypy:
                    error = (
                        "^Could not convert first line of 'file\.json' to "
                        "JSON \(Error when decoding null at char 1\)\. Line "
                        "is 'not JSON'\.$")
                else:
                    error = (
                        "^Could not convert first line of 'file\.json' to "
                        "JSON \(No JSON object could be decoded\)\. Line is "
                        "'not JSON'\.$")
            six.assertRaisesRegex(self, ValueError, error,
                                  DiamondReadsAlignments, reads, 'file.json',
                                  'database.fasta')

    def testScoreTitle_Bits(self):
        """
        The score title must be correct when we are using bit scores.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            self.assertEqual('Bit score', readsAlignments.params.scoreTitle)

    def testScoreTitle_EValue(self):
        """
        The score title must be correct when we are using e values.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = DiamondReadsAlignments(
                reads, 'file.json', 'database.fasta',
                scoreClass=LowerIsBetterScore)
            self.assertEqual('$- log_{10}(e)$',
                             readsAlignments.params.scoreTitle)

    def testNucleotides(self):
        """
        The nucleotide type of the subject must be correct.
        """
        params = deepcopy(PARAMS)
        mockOpener = mockOpen(read_data=dumps(params) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            self.assertFalse(readsAlignments.params.subjectIsNucleotides)

    def testApplicationParams(self):
        """
        DIAMOND parameters must be extracted from the input JSON file and
        stored correctly.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            self.assertEqual(PARAMS, readsAlignments.params.applicationParams)

    def testJSONParamsButNoHits(self):
        """
        When DIAMOND parameters are present in the input but there are no
        records, the __iter__ method of a L{DiamondReadsAlignments} instance
        must not yield anything.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            self.assertEqual([], list(readsAlignments))

    def testNotEnoughReads(self):
        """
        If a JSON file contains a parameters section and one hit, but there
        is no read to go with the hit, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(
            read_data=dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            error = ("Read generator failed to yield a read with id 'id0' as "
                     "found in record number 1 during parsing of DIAMOND "
                     "output file 'file\.json'\.")
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            six.assertRaisesRegex(self, ValueError, error, list,
                                  readsAlignments)

    def testIncorrectReadId(self):
        """
        If the query id of a hit does not match the id of the corresponding
        input read, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(
            read_data=dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('not id0', 'A' * 70))
            error = ("Read generator failed to yield a read with id 'id0' as "
                     "found in record number 1 during parsing of DIAMOND "
                     "output file 'file.json'.")
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            six.assertRaisesRegex(self, ValueError, error, list,
                                  readsAlignments)

    def testMoreReadsThanRecords(self):
        """
        If a JSON file contains a parameters section and one hit, but there
        are two query reads, the second read must still be returned, but have
        no alignments (i.e., length zero).
        """
        mockOpener = mockOpen(
            read_data=dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'G' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments)
            self.assertEqual(2, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)
            self.assertEqual(0, len(result[1]))

    def testOneJSONInput(self):
        """
        If a JSON file contains a parameters section and one record, it must
        be read correctly.
        """
        result = File([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])

        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.return_value = result
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            self.assertEqual(1, len(list(readsAlignments)))

    def testTwoJSONInputs(self):
        """
        If two JSON files are passed to L{DiamondReadsAlignments} each with a
        parameters section and one record, both records must be read correctly
        and the result should have 2 records.
        """

        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename, **kwargs):
                if self.first:
                    self.first = False
                    return File([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])
                else:
                    return File([dumps(PARAMS) + '\n', dumps(RECORD1) + '\n'])

        sideEffect = SideEffect()
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(
                reads, ['file1.json', 'file2.json'], 'database.fasta')
            result = list(readsAlignments)
            self.assertEqual(2, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)

    def testThreeJSONInputs(self):
        """
        If three JSON files are passed to L{DiamondReadsAlignments} with names
        that have a numeric prefix and each with a parameters section and one
        record, all records must be read correctly and the result should have
        3 records in the correct order.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('1.json', filename)
                    self.count += 1
                    return File([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])
                elif self.count == 1:
                    self.test.assertEqual('2.json', filename)
                    self.count += 1
                    return File([dumps(PARAMS) + '\n', dumps(RECORD1) + '\n'])
                else:
                    self.test.assertEqual('3.json', filename)
                    return File([dumps(PARAMS) + '\n', dumps(RECORD2) + '\n'])

        sideEffect = SideEffect(self)
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))

            # Note the files are given out of order. Their names will be
            # sorted before they are opened. The sorting of the names is
            # verified in the SideEffect class, above.
            readsAlignments = DiamondReadsAlignments(
                reads, ['3.json', '1.json', '2.json'], 'database.fasta')
            result = list(readsAlignments)
            self.assertEqual(3, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)
            self.assertEqual('id2', result[2].read.id)

    def testGetSubjectSequence(self):
        """
        The getSubjectSequence function must return an AAReadWithX instance
        with a string sequence.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, mode='r'):
                if self.count == 0:
                    self.test.assertEqual('file.json', filename)
                    self.count += 1
                    return File([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])
                elif self.count == 1:
                    self.count += 1
                    return File(['>id1 Description', 'AA\n'])
                else:
                    self.fail('Unexpected third call to open.')

        sideEffect = SideEffect(self)

        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = Reads()
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            subject = readsAlignments.getSubjectSequence('id1 Description')
            self.assertIsInstance(subject, AAReadWithX)
            self.assertIsInstance(subject.sequence, str)
            self.assertEqual('id1 Description', subject.id)
            self.assertEqual('AA', subject.sequence)

    def testHsps(self):
        """
        The hsps function must yield the HSPs.
        """
        # adjustHspsForPlotting changes HSPs in place, so we pass copied
        # records so we don't mess up other tests.
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n' +
            dumps(RECORD3) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            self.assertEqual(
                sorted([HSP(20), HSP(25), HSP(20), HSP(20), HSP(20), HSP(20)]),
                sorted(readsAlignments.hsps()))

    def testAdjustHspsForPlotting_EValueNoZero(self):
        """
        The adjustHspsForPlotting function must alter HSPs so that non-zero
        evalues are converted to the positive value of their negative exponent.
        """
        result = lambda a: File([
            dumps(PARAMS) + '\n', dumps(deepcopy(RECORD0)) + '\n',
            dumps(deepcopy(RECORD1)) + '\n', dumps(deepcopy(RECORD2)) + '\n',
            dumps(deepcopy(RECORD3)) + '\n'])

        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = result
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(
                reads, 'file.json', 'database.fasta',
                scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            title = 'gi|887699|gb|DQ37780 Cowpox virus 15'
            titleAlignments = titlesAlignments[title]
            readsAlignments.adjustHspsForPlotting(titleAlignments)
            hsps = sorted(titleAlignments.hsps())
            self.assertEqual([6.0, 5.0], [hsp.score.score for hsp in hsps])

    def testAdjustHspsForPlotting_EValueWithZero(self):
        """
        The adjustHspsForPlotting function must alter HSPs so that zero
        evalues are set randomly high.
        """
        result = lambda a: File([
            dumps(PARAMS) + '\n', dumps(deepcopy(RECORD0)) + '\n',
            dumps(deepcopy(RECORD1)) + '\n', dumps(deepcopy(RECORD2)) + '\n',
            dumps(deepcopy(RECORD3)) + '\n', dumps(deepcopy(RECORD4)) + '\n'])

        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = result
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            reads.add(Read('id3', 'A' * 70))
            reads.add(Read('id4', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(
                reads, 'file.json', 'database.fasta',
                scoreClass=LowerIsBetterScore)
            titlesAlignments = TitlesAlignments(readsAlignments)
            title = 'gi|887699|gb|DQ37780 Cowpox virus 15'
            titleAlignments = titlesAlignments[title]
            readsAlignments.adjustHspsForPlotting(titleAlignments)
            hsps = sorted(titleAlignments.hsps())
            # All we really know is that the first HSP will have a randomly
            # high value whose bounds we can check. The other values are
            # predictable.
            self.assertTrue(LSP(6.0 + 2) > hsps[0] >
                            LSP(6.0 + 2 + ZERO_EVALUE_UPPER_RANDOM_INCREMENT))
            self.assertEqual([6.0, 5.0, 3.0, 2.0],
                             [hsp.score.score for hsp in hsps[1:]])


class TestDiamondReadsAlignmentsFiltering(TestCase):
    """
    Test the DiamondReadsAlignments class filter function.
    """

    def testNoResultNoFilteringArgs(self):
        """
        If the L{DiamondReadsAlignments} filter function is called with no
        arguments, and there are no hits, it should produce a generator
        that yields no result.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter())
            self.assertEqual(0, len(result))

    def testOneHitNoFilteringArgs(self):
        """
        If the L{DiamondReadsAlignments} filter function is called with no
        arguments, and there is one hit, it should produce a generator that
        yields that hit.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter())
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)

    def testLimitZero(self):
        """
        If L{DiamondReadsAlignments} is limited to zero result, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(limit=0))
            self.assertEqual(0, len(result))

    def testLimitOne(self):
        """
        If L{DiamondReadsAlignments} is limited to one hit, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n' +
                              dumps(RECORD1) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(limit=1))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)

    def testOneAlignmentPerRead(self):
        """
        If L{DiamondReadsAlignments} is asked to deliver only the best
        alignment for each read, that must be respected.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 961,
                    "hsps": [
                        {
                            "sbjct_end": 869,
                            "expect": 1.25854e-10,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 836,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 150,
                            "btop": "",
                            "query_start": 362
                        },
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 740,
                    "hsps": [
                        {
                            "sbjct_end": 647,
                            "expect": 1.25e-43,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 614,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 180,
                            "btop": "",
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(oneAlignmentPerRead=True))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual('Merkel2', result[0][0].subjectTitle)

    def testScoreCutoffRemovesEntireAlignment_Bits(self):
        """
        If the L{DiamondReadsAlignments} filter function is supposed to filter
        on a scoreCutoff (bit score) and the cut-off value results in an
        alignment with no HSPs, then the alignment must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 961,
                    "hsps": [
                        {
                            "sbjct_end": 869,
                            "expect": 1.25854e-10,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 836,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 150,
                            "btop": "",
                            "query_start": 362
                        },
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 740,
                    "hsps": [
                        {
                            "sbjct_end": 647,
                            "expect": 1.25e-43,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 614,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 180,
                            "btop": "",
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(scoreCutoff=160))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual('Merkel2', result[0][0].subjectTitle)

    def testScoreCutoffRemovesEntireAlignment_EValue(self):
        """
        If the L{DiamondReadsAlignments} filter function is supposed to filter
        on a scoreCutoff (bit score) and the cut-off value results in an
        alignment with no HSPs, then the alignment must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 961,
                    "hsps": [
                        {
                            "sbjct_end": 869,
                            "expect": 1.25854e-10,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 836,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 150,
                            "btop": "",
                            "query_start": 362
                        },
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 740,
                    "hsps": [
                        {
                            "sbjct_end": 647,
                            "expect": 1.25e-43,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 614,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 180,
                            "btop": "",
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = DiamondReadsAlignments(
                reads, 'file.json', 'database.fasta',
                scoreClass=LowerIsBetterScore)
            result = list(readsAlignments.filter(scoreCutoff=1e-20))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual('Merkel2', result[0][0].subjectTitle)

    def testScoreCutoffRemovesHsps_Bits(self):
        """
        If the L{DiamondRecords} records function is supposed to filter on
        scoreCutoff (bit score) and the cut-off value results in some HSPs
        being invalid, then those HSPs must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 961,
                    "hsps": [
                        {
                            "sbjct_end": 869,
                            "expect": 1.25854e-10,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 836,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 150,
                            "btop": "",
                            "query_start": 362
                        },
                        {
                            "sbjct_end": 869,
                            "expect": 1.25e-20,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 836,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 170,
                            "btop": "",
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 740,
                    "hsps": [
                        {
                            "sbjct_end": 647,
                            "expect": 1.25e-43,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 614,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 180,
                            "btop": "",
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(scoreCutoff=160))

            # There should only be one HSP left in the alignments for the
            # first read, and it should have the right score (bit score).
            self.assertEqual(1, len(result[0][0].hsps))
            self.assertEqual(HSP(170), result[0][0].hsps[0])

            # The second alignment should also be present.
            self.assertEqual(1, len(result[0][1].hsps))
            self.assertEqual(HSP(180), result[0][1].hsps[0])

    def testScoreCutoffRemovesHsps_EValue(self):
        """
        If the L{DiamondRecords} records function is supposed to filter on
        scoreCutoff (bit score) and the cut-off value results in some HSPs
        being invalid, then those HSPs must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 961,
                    "hsps": [
                        {
                            "sbjct_end": 869,
                            "expect": 1.25854e-10,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 836,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 150,
                            "btop": "",
                            "query_start": 362
                        },
                        {
                            "sbjct_end": 869,
                            "expect": 1.25e-20,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 836,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 170,
                            "btop": "",
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 740,
                    "hsps": [
                        {
                            "sbjct_end": 647,
                            "expect": 1.25e-43,
                            "sbjct": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "sbjct_start": 614,
                            "query": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                      "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
                            "frame": 1,
                            "query_end": 462,
                            "bits": 180,
                            "btop": "",
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = DiamondReadsAlignments(
                reads, 'file.json', 'database.fasta',
                scoreClass=LowerIsBetterScore)
            result = list(readsAlignments.filter(scoreCutoff=1e-15))

            # There should only be one HSP left in the alignments for the
            # first read, and it should have the right score (e-value).
            self.assertEqual(1, len(result[0][0].hsps))
            self.assertEqual(LSP(1.25e-20), result[0][0].hsps[0])

            # The second alignment should also be present.
            self.assertEqual(1, len(result[0][1].hsps))
            self.assertEqual(LSP(1.25e-43), result[0][1].hsps[0])

    def testTitleByRegexCaseInvariant(self):
        """
        Filtering with a title regex must work independent of case.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(titleRegex='sqUIRRel'))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0][0].subjectTitle)

    def testTitleByRegexAllAlignments(self):
        """
        Filtering with a title regex must work in the case that all alignments
        for a hit match the regex.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(titleRegex='squirrel'))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0][0].subjectTitle)

    def testTitleByRegexOneAlignments(self):
        """
        Filtering with a title regex must work in the case that only some
        alignments for a hit match the regex.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(titleRegex='Mummy'))
            self.assertEqual(1, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual('gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.',
                             result[0][0].subjectTitle)

    def testTitleByNegativeRegexOneAlignment(self):
        """
        Filtering with a negative title regex must work in the case that only
        some alignments for a hit are ruled out (in which case only those
        alignments must be removed but the hit is still valid).
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(negativeTitleRegex='Mummy'))
            self.assertEqual(3, len(result))
            self.assertEqual('id1', result[1].read.id)
            self.assertEqual(1, len(result[1]))
            self.assertEqual('gi|887699|gb|DQ37780 Monkeypox virus 456',
                             result[1][0].subjectTitle)

    def testTitleByNegativeRegexMatchesAll(self):
        """
        Filtering with a negative title regex that matches all alignments
        must remove everything and return an empty result.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(negativeTitleRegex='pox'))
            self.assertEqual(0, len(result))

    def testTitleByNegativeRegexMatchingAllWithWhitelist(self):
        """
        Filtering with a negative title regex that matches all alignments
        must remove everything and result in no hits, except for any
        whitelisted titles.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            title = 'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99'
            result = list(readsAlignments.filter(negativeTitleRegex='pox',
                                                 whitelist=[title]))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual(title, result[0][0].subjectTitle)

    def testTitleByRegexMatchingAllWithBlacklist(self):
        """
        Filtering with a title regex that matches all alignments
        must keep everything, except for any blacklisted titles.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            blacklist = ['gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                         'gi|887699|gb|DQ37780 Squirrelpox virus 55']
            result = list(readsAlignments.filter(titleRegex='pox',
                                                 blacklist=blacklist))
            self.assertEqual(2, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual('id2', result[1].read.id)

    def testTitleTruncation(self):
        """
        When truncating titles, if a set of matched sequences has titles that
        are identical up to the truncation word, only the first found is
        returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = readsAlignments.filter(truncateTitlesAfter='virus')
            result = list(result)
            self.assertEqual(3, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            # The Squirrelpox virus 55 hit in RECORD0 is not returned.
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0][0].subjectTitle)

    def testMinTitleSequenceLength(self):
        """
        It must be possible to filter alignments based on minimum hit sequence
        length.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(minSequenceLen=37500))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 55',
                             result[0][0].subjectTitle)

    def testMinTitleSequenceLengthNoHits(self):
        """
        It must be possible to filter alignments based on minimum hit sequence
        length and if nothing sufficiently long matches, an empty list of
        alignments must be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(minSequenceLen=1000000))
            self.assertEqual(0, len(result))

    def testMaxTitleSequenceLength(self):
        """
        It must be possible to filter alignments based on maximum hit sequence
        length.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(maxSequenceLen=31000))
            self.assertEqual(1, len(result))
            self.assertEqual('id2', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual('gi|887699|gb|DQ37780 Cowpox virus 15',
                             result[0][0].subjectTitle)

    def testMaxTitleSequenceLengthNoHits(self):
        """
        It must be possible to filter alignments based on maximum hit sequence
        length and if no sufficiently short sequences match, an empty
        list of alignments must be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(maxSequenceLen=10000))
            self.assertEqual(0, len(result))

    def testMinAndMaxTitleSequenceLength(self):
        """
        It must be possible to filter alignments simultaneously on minimum and
        maximum hit sequence length.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(minSequenceLen=37000,
                                                 maxSequenceLen=38000))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(2, len(result[0]))
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0][0].subjectTitle)
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 55',
                             result[0][1].subjectTitle)

    def testMinStart(self):
        """
        It must be possible to filter alignments based on minimum offset in
        the hit sequence.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(minStart=15300))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual('gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                             result[0][0].subjectTitle)

    def testMinStartNoHits(self):
        """
        It must be possible to filter alignments based on minimum offset in
        the hit sequence, and if no hsps match then an empty result set
        must be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(minStart=100000))
            self.assertEqual(0, len(result))

    def testMaxStop(self):
        """
        It must be possible to filter alignments based on maximum offset in
        the hit sequence.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(maxStop=1500))
            self.assertEqual(1, len(result))
            self.assertEqual('id2', result[0].read.id)
            self.assertEqual(1, len(result[0]))
            self.assertEqual('gi|887699|gb|DQ37780 Cowpox virus 15',
                             result[0][0].subjectTitle)

    def testMaxStopNoHits(self):
        """
        It must be possible to filter alignments based on maximum offset in
        the hit sequence, and if no hsps match then an empty result set must
        be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(maxStop=100))
            self.assertEqual(0, len(result))

    def testMinStartAndMaxstop(self):
        """
        It must be possible to filter alignments based simultaneously on
        mininum and maximum offset in the hit sequence.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(minStart=9000, maxStop=12000))
            self.assertEqual(1, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual(2, len(result[0]))

    def testRepeatedFilter_MinStartThenMinStart(self):
        """
        It must be possible to filter alignments multiple times using the same
        filter parameters.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            readsAlignments.filter(minStart=9000)
            readsAlignments.filter(minStart=9000)
            result = list(readsAlignments)
            self.assertEqual(2, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)

    def testRepeatedFilter_MinStartThenMaxstop(self):
        """
        It must be possible to filter alignments multiple times using different
        filter parameters.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            readsAlignments.filter(minStart=9000)
            readsAlignments.filter(maxStop=12000)
            result = list(readsAlignments)
            self.assertEqual(1, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual(2, len(result[0]))

    def testClearFilter(self):
        """
        It must be possible to clear any filtering that has been applied.
        """
        fp = File([
            dumps(PARAMS) + '\n', dumps(RECORD0) + '\n',
            dumps(RECORD1) + '\n', dumps(RECORD2) + '\n'])
        with patch.object(builtins, 'open') as mockOpen:
            mockOpen.return_value = fp
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            self.assertEqual(3, len(list(readsAlignments)))
            readsAlignments.filter(minStart=9000)
            readsAlignments.filter(maxStop=12000)
            readsAlignments.filter(scoreCutoff=19)
            result = list(readsAlignments)
            self.assertEqual(1, len(result))
            readsAlignments.clearFilter()
            self.assertEqual(3, len(list(readsAlignments)))

    def testReadIdNoMatches(self):
        """
        When filtering on alignments based on a regex for
        read ids that matches no ids, an empty generator must be returned.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(readIdRegex='blah'))
            self.assertEqual(0, len(result))

    def testReadId(self):
        """
        It must be possible to filter alignments based on a regex for
        read ids.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(readIdRegex='id[12]'))
            self.assertEqual(2, len(result))
            self.assertEqual('id1', result[0].read.id)
            self.assertEqual('id2', result[1].read.id)

    def testReadIdAnchored(self):
        """
        It must be possible to filter alignments based on a regex for
        read ids that is anchored at start and end.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(readIdRegex='^id0$'))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)

    def testReadIdCaseSensitive(self):
        """
        Filtering alignments based on a regex for read ids must be case
        sensitive.
        """
        mockOpener = mockOpen(read_data=(
            dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n' +
            dumps(RECORD1) + '\n' + dumps(RECORD2) + '\n'))
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = DiamondReadsAlignments(reads, 'file.json',
                                                     'database.fasta')
            result = list(readsAlignments.filter(readIdRegex='^ID0$'))
            self.assertEqual(0, len(result))

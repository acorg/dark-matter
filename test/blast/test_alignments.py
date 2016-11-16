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

from six import StringIO
from Bio import SeqIO

from ..mocking import mockOpen, File
from .sample_data import PARAMS, RECORD0, RECORD1, RECORD2, RECORD3, RECORD4

from dark.reads import Read, Reads, DNARead
from dark.hsp import HSP, LSP
from dark.score import LowerIsBetterScore
from dark.blast.alignments import (
    BlastReadsAlignments,  ZERO_EVALUE_UPPER_RANDOM_INCREMENT)
from dark.titles import TitlesAlignments
from dark import ncbidb


class TestBlastReadsAlignments(TestCase):
    """
    Test the BlastReadsAlignments class.
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
                                  BlastReadsAlignments, reads, 'file.json')

    def testNonJSONInput(self):
        """
        When given a file whose contents are not JSON, attempting to
        read the BLAST hits from it must raise a C{ValueError}.
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
                                  BlastReadsAlignments, reads, 'file.json')

    def testScoreTitle_Bits(self):
        """
        The score title must be correct when we are using bit scores.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual('Bit score', readsAlignments.params.scoreTitle)

    def testScoreTitle_EValue(self):
        """
        The score title must be correct when we are using e values.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            self.assertEqual('$- log_{10}(e)$',
                             readsAlignments.params.scoreTitle)

    def testNucleotidesBlastn(self):
        """
        The nucleotide type of the subject must be correct when we are using
        blastn.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual('blastn', readsAlignments.params.application)
            self.assertTrue(readsAlignments.params.subjectIsNucleotides)

    def testNucleotidesTblastx(self):
        """
        The nucleotide type of the subject must be correct when we are using
        tblastx.
        """
        params = deepcopy(PARAMS)
        params['application'] = 'tblastx'
        mockOpener = mockOpen(read_data=dumps(params) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual('tblastx', readsAlignments.params.application)
            self.assertTrue(readsAlignments.params.subjectIsNucleotides)

    def testNucleotidesBlastx(self):
        """
        The nucleotide type of the subject must be correct when we are using
        blastx.
        """
        params = deepcopy(PARAMS)
        params['application'] = 'blastx'
        mockOpener = mockOpen(read_data=dumps(params) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual('blastx', readsAlignments.params.application)
            self.assertFalse(readsAlignments.params.subjectIsNucleotides)

    def testApplicationParams(self):
        """
        BLAST parameters must be extracted from the input JSON file and stored
        correctly.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual(PARAMS, readsAlignments.params.applicationParams)

    def testJSONParamsButNoHits(self):
        """
        When BLAST parameters are present in the input but there are no
        records, the __iter__ method of a L{BlastReadsAlignments} instance must
        not yield anything.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            error = ("Read generator failed to yield read number 1 during "
                     "parsing of BLAST file 'file\.json'\.")
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            six.assertRaisesRegex(self, ValueError, error, list,
                                  readsAlignments)

    def testTooManyReads(self):
        """
        If a JSON file contains a parameters section and one hit, but there
        is more than one read, a C{ValueError} must be raised.
        """
        mockOpener = mockOpen(
            read_data=dumps(PARAMS) + '\n' + dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'G' * 70))
            error = ("Reads iterator contained more reads than the number of "
                     "BLAST records found \(1\)\. First unknown read id is "
                     "'id1'\.")
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            error = ("The reads you have provided do not match the BLAST "
                     "output: BLAST record query id \(id0\) does "
                     "not match the id of the supposedly corresponding read "
                     "\(not id0\)\.")
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            six.assertRaisesRegex(self, ValueError, error, list,
                                  readsAlignments)

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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual(1, len(list(readsAlignments)))

    def testTwoJSONInputs(self):
        """
        If two JSON files are passed to L{BlastReadsAlignments} each with a
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
            readsAlignments = BlastReadsAlignments(
                reads, ['file1.json', 'file2.json'])
            result = list(readsAlignments)
            self.assertEqual(2, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)

    def testThreeJSONInputs(self):
        """
        If three JSON files are passed to L{BlastReadsAlignments} with names
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
            readsAlignments = BlastReadsAlignments(
                reads, ['3.json', '1.json', '2.json'])
            result = list(readsAlignments)
            self.assertEqual(3, len(result))
            self.assertEqual('id0', result[0].read.id)
            self.assertEqual('id1', result[1].read.id)
            self.assertEqual('id2', result[2].read.id)

    def testIncompatibleParameters(self):
        """
        If two compressed (bz2) JSON files with incompatible parameters
        are given to L{BlastReadsAlignments}, a C{ValueError} must be
        raised when the files are read.
        """

        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return File([dumps(PARAMS) + '\n', dumps(RECORD0) + '\n'])
                else:
                    params = deepcopy(PARAMS)
                    params['application'] = 'Skype'
                    return File([dumps(params) + '\n', dumps(RECORD1) + '\n'])

        sideEffect = SideEffect()
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            if six.PY3:
                error = (
                    "^Incompatible BLAST parameters found\. The parameters "
                    "in file2\.json differ from those originally found "
                    "in file1\.json. Summary of differences:\n\tParam "
                    "'application' initial value 'BLASTN' differs from "
                    "later value 'Skype'$")
            else:
                # Python 2 prints a 'u' before the repr of strings in the error
                # message. In Python 3 all strings are unicode.
                error = (
                    "^Incompatible BLAST parameters found\. The parameters "
                    "in file2\.json differ from those originally found "
                    "in file1\.json. Summary of differences:\n\tParam "
                    "u'application' initial value u'BLASTN' differs from "
                    "later value u'Skype'$")
            readsAlignments = BlastReadsAlignments(
                reads, ['file1.json', 'file2.json'])
            six.assertRaisesRegex(self, ValueError, error, list,
                                  readsAlignments)

    def testGetSubjectSequence(self):
        """
        The getSubjectSequence function must return a correct C{DNARead}
        instance with a 'sequence' attribute that is a string.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            with patch.object(ncbidb, 'getSequence') as mockMethod:
                mockMethod.return_value = SeqIO.read(
                    StringIO('>id1 Description\nAA\n'), 'fasta')
                sequence = readsAlignments.getSubjectSequence('title')
                self.assertIsInstance(sequence, DNARead)
                self.assertIsInstance(sequence.sequence, str)
                self.assertEqual('id1 Description', sequence.id)
                self.assertEqual('AA', sequence.sequence)

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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            self.assertEqual(
                sorted([HSP(20), HSP(25), HSP(20), HSP(20), HSP(20), HSP(20)]),
                sorted(readsAlignments.hsps()))

    def testAdjustHspsForPlotting_EValueNoZero(self):
        """
        The adjustHspsForPlotting function must alter HSPs so that non-zero
        evalues are converted to the positive value of their negative exponent.
        """
        result = lambda a, **kwargs: File([
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
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
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
        result = lambda a, **kwargs: File([
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
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
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


class TestBlastReadsAlignmentsFiltering(TestCase):
    """
    Test the BlastReadsAlignments class filter function.
    """

    def testNoResultNoFilteringArgs(self):
        """
        If the L{BlastReadsAlignments} filter function is called with no
        arguments, and there are no hits, it should produce a generator
        that yields no result.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter())
            self.assertEqual(0, len(result))

    def testOneHitNoFilteringArgs(self):
        """
        If the L{BlastReadsAlignments} filter function is called with no
        arguments, and there is one hit, it should produce a generator that
        yields that hit.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter())
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)

    def testLimitZero(self):
        """
        If L{BlastReadsAlignments} is limited to zero result, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(limit=0))
            self.assertEqual(0, len(result))

    def testLimitOne(self):
        """
        If L{BlastReadsAlignments} is limited to one hit, that limit must
        be respected.
        """
        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(RECORD0) + '\n' +
                              dumps(RECORD1) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(limit=1))
            self.assertEqual(1, len(result))
            self.assertEqual('id0', result[0].read.id)

    def testOneAlignmentPerRead(self):
        """
        If L{BlastReadsAlignments} is asked to deliver only the best alignment
        for each read, that must be respected.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25854e-43,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 182.092,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 2220,
                    "hsps": [
                        {
                            "sbjct_end": 1841,
                            "expect": 1.25854e-43,
                            "sbjct": "AATCCAGGGAATCTAATAAAATAATCAA",
                            "sbjct_start": 1942,
                            "query": "AATCCAGGGAATCTTAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 180,
                            "query_start": 362
                        }
                    ],
                    "title":"Merkel2"
                }
            ]
        }

        mockOpener = mockOpen(read_data=dumps(PARAMS) + '\n' +
                              dumps(record) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            reads = Reads()
            reads.add(Read('H6E8I1T01BFUH9', 'A' * 500))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(oneAlignmentPerRead=True))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual('Merkel1', result[0][0].subjectTitle)

    def testScoreCutoffRemovesEntireAlignment_Bits(self):
        """
        If the L{BlastReadsAlignments} filter function is supposed to filter on
        a scoreCutoff (bit score) and the cut-off value results in an
        alignment with no HSPs, then the alignment must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25854e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 150,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 2220,
                    "hsps": [
                        {
                            "sbjct_end": 1841,
                            "expect": 1.25854e-43,
                            "sbjct": "AATCCAGGGAATCTAATAAAATAATCAA",
                            "sbjct_start": 1942,
                            "query": "AATCCAGGGAATCTTAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 180,
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(scoreCutoff=160))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual('Merkel2', result[0][0].subjectTitle)

    def testScoreCutoffRemovesEntireAlignment_EValue(self):
        """
        If the L{BlastReadsAlignments} filter function is supposed to filter on
        a scoreCutoff (bit score) and the cut-off value results in an
        alignment with no HSPs, then the alignment must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 150,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 2220,
                    "hsps": [
                        {
                            "sbjct_end": 1841,
                            "expect": 1.25e-30,
                            "sbjct": "AATCCAGGGAATCTAATAAAATAATCAA",
                            "sbjct_start": 1942,
                            "query": "AATCCAGGGAATCTTAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 180,
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
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            result = list(readsAlignments.filter(scoreCutoff=1e-20))
            self.assertEqual(1, len(result))
            self.assertEqual(1, len(result[0]))
            self.assertEqual('Merkel2', result[0][0].subjectTitle)

    def testScoreCutoffRemovesHsps_Bits(self):
        """
        If the L{BlastRecords} records function is supposed to filter on
        scoreCutoff (bit score) and the cut-off value results in some HSPs
        being invalid, then those HSPs must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25854e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 150,
                            "query_start": 362
                        },
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25e-20,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 170,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 2220,
                    "hsps": [
                        {
                            "sbjct_end": 1841,
                            "expect": 1.25e-43,
                            "sbjct": "AATCCAGGGAATCTAATAAAATAATCAA",
                            "sbjct_start": 1942,
                            "query": "AATCCAGGGAATCTTAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 180,
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
        If the L{BlastRecords} records function is supposed to filter on
        scoreCutoff (bit score) and the cut-off value results in some HSPs
        being invalid, then those HSPs must be removed entirely.
        """
        record = {
            "query": "H6E8I1T01BFUH9",
            "alignments": [
                {
                    "length": 2885,
                    "hsps": [
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25e-10,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 150,
                            "query_start": 362
                        },
                        {
                            "sbjct_end": 2506,
                            "expect": 1.25e-20,
                            "sbjct": "AATCCAGGGAATGAATAAAATAATCATTAGCAGTAACAA",
                            "sbjct_start": 2607,
                            "query": "AATCCAGGGAATAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 170,
                            "query_start": 362
                        }
                    ],
                    "title": "Merkel1"
                },
                {
                    "length": 2220,
                    "hsps": [
                        {
                            "sbjct_end": 1841,
                            "expect": 1.25e-30,
                            "sbjct": "AATCCAGGGAATCTAATAAAATAATCAA",
                            "sbjct_start": 1942,
                            "query": "AATCCAGGGAATCTTAAA-TAATCATTAGCAGTAACAA",
                            "frame": [1, -1],
                            "query_end": 462,
                            "bits": 180,
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
            readsAlignments = BlastReadsAlignments(
                reads, 'file.json', scoreClass=LowerIsBetterScore)
            result = list(readsAlignments.filter(scoreCutoff=1e-15))

            # There should only be one HSP left in the alignments for the
            # first read, and it should have the right score (e-value).
            self.assertEqual(1, len(result[0][0].hsps))
            self.assertEqual(LSP(1.25e-20), result[0][0].hsps[0])

            # The second alignment should also be present.
            self.assertEqual(1, len(result[0][1].hsps))
            self.assertEqual(LSP(1.25e-30), result[0][1].hsps[0])

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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
        result = lambda a: File([
            dumps(PARAMS) + '\n', dumps(RECORD0) + '\n',
            dumps(RECORD1) + '\n', dumps(RECORD2) + '\n'])

        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = result
            reads = Reads()
            reads.add(Read('id0', 'A' * 70))
            reads.add(Read('id1', 'A' * 70))
            reads.add(Read('id2', 'A' * 70))
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
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
            readsAlignments = BlastReadsAlignments(reads, 'file.json')
            result = list(readsAlignments.filter(readIdRegex='^ID0$'))
            self.assertEqual(0, len(result))

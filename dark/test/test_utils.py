from unittest import TestCase
from dark.utils import findCodons, normalizeHSP, interestingRecords
from Bio.Seq import Seq


class HSP(object):
    def __init__(self, subjectStart, subjectEnd, queryStart, queryEnd):
        """
        A fake HSP class (complete with 1-based offsets).
        """
        self.sbjct_start = subjectStart
        self.sbjct_end = subjectEnd
        self.query_start = queryStart
        self.query_end = queryEnd
        self.frame = (1, 1 if subjectStart < subjectEnd else -1)


class NormalizeHSPSubjectOffsetsAscending(TestCase):
    """
    Tests for the utils.normalizeHSP function when the subject start is less
    than the subject end.
    """

    def testSubjectAndQueryIdentical(self):
        """The subject start and end are identical to those of the query.

             ssss
             qqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=4)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': 0,
            'queryEnd': 4,
        }, normalized)

    def testSubjectOverlapsQueryToTheRight(self):
        """The subject overlaps the query to the right.

              ssssss
              qqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=6, queryStart=1, queryEnd=4)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 6,
            'queryStart': 0,
            'queryEnd': 4,
        }, normalized)

    def testSubjectOverlapsQueryToTheLeft(self):
        """The subject overlaps the query to the left.

              ssssss
                qqqq
        """
        hsp = HSP(subjectStart=3, subjectEnd=6, queryStart=1, queryEnd=4)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'queryStart': 2,
            'queryEnd': 6,
        }, normalized)

    def testQueryExtendsLeft(self):
        """
        The query sticks out to the left of the subject.

               ssss
             qqqqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=6)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': -2,
            'queryEnd': 4,
        }, normalized)

    def testQueryExtendsRight(self):
        """The query sticks out to the right of the subject.

               ssss
               qqqqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=6)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': 0,
            'queryEnd': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """The query extends to the right and left of the subject.

                ssss
               qqqqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=5)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': -1,
            'queryEnd': 5,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """The subject extends to the right and left of the query.

               sssssss
                qqqq
        """
        hsp = HSP(subjectStart=2, subjectEnd=5, queryStart=1, queryEnd=4)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'queryStart': 1,
            'queryEnd': 5,
        }, normalized)


class NormalizeHSPSubjectOffsetsDescending(TestCase):
    """
    Tests for the utils.normalizeHSP function when the subject start is greater
    than the subject end.
    """

    def testSubjectAndQueryIdentical(self):
        """
        The subject start and end are identical to those of the query.

              ssss
              qqqq
        """
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=1, queryEnd=4)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': 0,
            'queryEnd': 4,
        }, normalized)

    def testQueryExtendsLeft(self):
        """The query sticks out to the left of the subject.

                ssss
                qqqqqq
        """
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=1, queryEnd=4)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': -2,
            'queryEnd': 4,
        }, normalized)

    def testQueryExtendsLeft2(self):
        """The query sticks out to the left of the subject.

              ssss
                qqqqqq
        """
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=3, queryEnd=6)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': 0,
            'queryEnd': 6,
        }, normalized)

    def testQueryExtendsRight(self):
        """The query sticks out to the right of the subject.

                  ssss
                qqqqqq
        """
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=3, queryEnd=6)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': 0,
            'queryEnd': 6,
        }, normalized)

    def testQueryExtendsRight2(self):
        """The query sticks out to the right of the query.

                ssss
            qqqqqq
        """
        hsp = HSP(subjectStart=2, subjectEnd=1, queryStart=5, queryEnd=6)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 2,
            'queryStart': 0,
            'queryEnd': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """The query extends to the right and left of the subject.

                  ssss
                qqqqqqq
        """
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=3, queryEnd=6)
        normalized = normalizeHSP(hsp, 7)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': -1,
            'queryEnd': 6,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """The subject extends to the right and left of the query.

                sssssss
                 qqqq
        """
        hsp = HSP(subjectStart=5, subjectEnd=2, queryStart=1, queryEnd=4)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'queryStart': 1,
            'queryEnd': 5,
        }, normalized)


class TestFindCodons(TestCase):
    """Tests of the findCodons helper. """

    def testNoMatches(self):
        """
        When there are no codons in the sequence, returns an empty list.
        """
        seq = Seq('AAAAAA')
        self.assertEqual([], list(findCodons(seq, set(['ATG', 'AGG']))))

    def testMatchAtStart(self):
        """
        Finds a codon at the start of the sequence.
        """
        seq = Seq('ATGAAA')
        self.assertEqual([0], list(findCodons(seq, set(['ATG', 'AGG']))))

    def testMatchAtEnd(self):
        """
        Finds a codon at the end of the sequence.
        """
        seq = Seq('ATGAAA')
        self.assertEqual([3], list(findCodons(seq, set(['AAA']))))

    def testMatchMultiple(self):
        """
        Finds multiple codons in the sequence.
        """
        seq = Seq('ATGAAAGGGCCC')
        self.assertEqual([0, 9], list(findCodons(seq, set(['ATG', 'CCC']))))

    def testDoesNotFindOutOfFrameMatches(self):
        """
        Does not find matching codons that are in a non-zero frame in the
        sequence.
        """
        seq = Seq('TATGAAAGGGCCC')
        self.assertEqual([], list(findCodons(seq, set(['ATG', 'CCC']))))


class InterestingRecords(TestCase):
    """
    Tests for the utils.interestingRecords function.
    """

    def testNoSequences(self):
        """
        If we give interestingRecords no records, we should get nothing back.
        """
        self.assertEqual({}, interestingRecords({}))

    def testOneSequenceWithNoFiltering(self):
        """
        If we give interestingRecords one record, but no filtering restrictions
        we should get the same record back
        """
        self.assertEqual({'a': 23}, interestingRecords({'a': 23}))

    def testFilterTitle(self):
        """
        Filtering with a title (i.e., dictionary key) should work.
        """
        self.assertEqual(
            {
                'a': 23
            },
            interestingRecords(
                {
                    'a': 23,
                    'b': 15
                },
                keyRegex='a'
            ))

    def testFilterLengthMinMax(self):
        """
        Filtering with a minimum and maximum sequence length should work.
        """
        self.assertEqual(
            {
                'b': {
                    'length': 20
                }
            },
            interestingRecords(
                {
                    'a': {
                        'length': 10
                    },
                    'b': {
                        'length': 20
                    },
                    'c': {
                        'length': 30
                    }
                },
                minSequenceLen=15,
                maxSequenceLen=25
            ))

    def testFilterMinMatchingReads(self):
        """
        Filtering with a minimum number of reads should work.
        """
        self.assertEqual(
            {
                'b': {
                    'count': 20
                },
                'c': {
                    'count': 30
                }
            },
            interestingRecords(
                {
                    'a': {
                        'count': 10
                    },
                    'b': {
                        'count': 20
                    },
                    'c': {
                        'count': 30
                    }
                },
                minMatchingReads=15
            ))

    def testFilterMaxMeanEValue(self):
        """
        Filtering with a maximum mean e-value should work.
        """
        self.assertEqual(
            {
                'b': {
                    'eMean': 5
                }
            },
            interestingRecords(
                {
                    'a': {
                        'eMean': 10
                    },
                    'b': {
                        'eMean': 5
                    },
                    'c': {
                        'eMean': 30
                    }
                },
                maxMeanEValue=6
            ))

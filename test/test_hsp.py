from unittest import TestCase
from dark.hsp import normalizeHSP


class Frame(object):
    def __init__(self, query, subject):
        self.query = query
        self.subject = subject


class HSP(object):
    def __init__(self, subjectStart, subjectEnd, queryStart, queryEnd, frame,
                 subject='', query=''):
        """
        A fake HSP class (with 1-based offsets, as are used in BLAST).
        """
        self.sbjct_start = subjectStart
        self.sbjct_end = subjectEnd
        self.query_start = queryStart
        self.query_end = queryEnd
        self.frame = (frame.query, frame.subject)
        self.sbjct = subject
        self.query = query

        # The following assertion is not valid.
        # assert abs(subjectEnd - subjectStart) == abs(queryEnd - queryStart)
        # That's because BLAST might find a match that requires a gap in the
        # query or subject. The indices that it reports do not include the gap
        # and so the differences in the lengths of the sections of the query
        # and subject may not be the same.


class Template(object):
    def __init__(self, template):
        template = template.split('\n')
        # Allow the first template line to be empty.
        if len(template[0]) == 0:
            template = template[1:]
        self.subject = template[0].rstrip()
        self.query = template[1].rstrip()
        self.queryResult = template[2].rstrip()

        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen, positive) = \
            self.analyze(self.subject)
        origin = spacesLen

        self.matchLen = matchLen
        self.subjectStart = 0
        self.subjectLen = len(self.subject) - spacesLen
        self.subjectMatchStart = leadingDotsLen
        self.subjectPositive = positive

        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen, positive) = \
            self.analyze(self.query)

        assert self.matchLen == matchLen
        self.queryStart = spacesLen - origin
        self.queryLen = len(self.query) - spacesLen
        self.queryMatchStart = leadingDotsLen
        self.queryPositive = positive

        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen, positive) = \
            self.analyze(self.queryResult)

        assert self.matchLen == matchLen
        self.queryResultStart = spacesLen - origin
        self.queryResultLen = len(self.queryResult) - spacesLen
        assert self.queryResultLen == self.queryLen
        self.queryResultMatchStart = leadingDotsLen

    def leadingCharMatchLen(self, str, chars=' '):
        return len(str) - len(str.lstrip(chars))

    def analyze(self, str):
        offset = spacesLen = self.leadingCharMatchLen(str)
        leadingDotsLen = self.leadingCharMatchLen(str[offset:], '.')
        offset += leadingDotsLen
        assert str[offset] in ('>', '<'), 'Oops: "%s"' % str[leadingDotsLen]
        positive = str[offset] == '>'
        matchLen = self.leadingCharMatchLen(str[offset:], '<>')
        offset += matchLen
        trailingDotsLen = self.leadingCharMatchLen(str[offset:], '.')
        return (spacesLen, leadingDotsLen, matchLen, trailingDotsLen, positive)

    def hsp(self):
        """
        Make an HSP. Use 1-based offsets.
        """
        return HSP(subjectStart=self.subjectMatchStart + 1,
                   subjectEnd=self.subjectMatchStart + self.matchLen,
                   queryStart=self.queryMatchStart + 1,
                   queryEnd=self.queryMatchStart + self.matchLen,
                   frame=Frame(
                       1 if self.queryPositive else -1,
                       1 if self.subjectPositive else -1,
                   ))


class TestTemplate(TestCase):
    """
    Tests for our helper Template class.
    """

    def testIt(self):
        template = Template('''
                                             ....>>>...........
                                               ..<<<...............
                                  ...............>>>..
        ''')
        self.assertEqual(3, template.matchLen)

        self.assertTrue(template.subjectPositive)
        self.assertEqual(0, template.subjectStart)
        self.assertEqual(18, template.subjectLen)
        self.assertEqual(4, template.subjectMatchStart)

        self.assertFalse(template.queryPositive)
        self.assertEqual(2, template.queryStart)
        self.assertEqual(20, template.queryLen)
        self.assertEqual(2, template.queryMatchStart)

        self.assertEqual(20, template.queryResultLen)
        self.assertEqual(-11, template.queryResultStart)
        self.assertEqual(15, template.queryResultMatchStart)

        hsp = template.hsp()
        self.assertEqual(5, hsp.sbjct_start)
        self.assertEqual(7, hsp.sbjct_end)
        self.assertEqual(3, hsp.query_start)
        self.assertEqual(5, hsp.query_end)
        self.assertEqual((-1, 1), hsp.frame)


class TestNormalizeHSPMixin(object):
    """
    Tests for normalizeHSP when the query and subject are both positive.
    """

    def check(self, templateStr):
        template = Template(templateStr)
        normalized = normalizeHSP(template.hsp(), template.queryLen)
        self.assertEqual({
            'subjectStart': template.subjectMatchStart,
            'subjectEnd': template.subjectMatchStart + template.matchLen,
            'queryStart': template.queryResultStart,
            'queryEnd': template.queryResultStart + template.queryLen,
        }, normalized)


class SubjectPositiveQueryPositive(TestCase, TestNormalizeHSPMixin):
    def testIdentical1(self):
        self.check('''
                                    >
                                    >
                                    >
                  ''')

    def testIdentical2(self):
        self.check('''
                                    ....>>>>...
                                    ....>>>>...
                                    ....>>>>...
                  ''')

    def testSubjectExtendsLeft1(self):
        self.check('''
                                    ....>>
                                        >>
                                        >>
                  ''')

    def testSubjectExtendsLeft2(self):
        self.check('''
                                    ....>>>...
                                        >>>...
                                        >>>...
                  ''')

    def testSubjectExtendsRight1(self):
        self.check('''
                                        >>......
                                        >>...
                                        >>...
                  ''')

    def testSubjectExtendsRight2(self):
        self.check('''
                                        ..>>>>......
                                        ..>>>>...
                                        ..>>>>...
                  ''')

    def testSubjectExtendsBoth(self):
        self.check('''
                                    ....>>>...........
                                      ..>>>....
                                      ..>>>....
                  ''')

    def testQueryExtendsLeft1(self):
        self.check('''
                                        >>
                                    ....>>
                                    ....>>
                  ''')

    def testQueryExtendsLeft2(self):
        self.check('''
                                        >>>...
                                    ....>>>...
                                    ....>>>...
                  ''')

    def testQueryExtendsRight1(self):
        self.check('''
                                        >>...
                                        >>......
                                        >>......
                  ''')

    def testQueryExtendsRight2(self):
        self.check('''
                                        ..>>>>...
                                        ..>>>>......
                                        ..>>>>......
                  ''')

    def testQueryExtendsBoth(self):
        self.check('''
                                      ..>>>....
                                    ....>>>...........
                                    ....>>>...........
                  ''')

    def testSubjectExtendsLeftQueryExtendsRight(self):
        self.check('''
                                    ....>>>...........
                                      ..>>>...............
                                      ..>>>...............
                  ''')

    def testSubjectExtendsRightQueryExtendsLeft(self):
        self.check('''
                                      ..>>>...............
                                    ....>>>...........
                                    ....>>>...........
                  ''')


class SubjectPositiveQueryNegative(TestCase, TestNormalizeHSPMixin):

    def testIdentical1(self):
        self.check('''
                                    >
                                    <
                                    >
                  ''')

    def testIdentical2(self):
        self.check('''
                                    >>>>
                                    <<<<
                                    >>>>
                  '''),

    def testSubjectExtendsLeft1(self):
        self.check('''
                                    ....>>>
                                      ..<<<
                                        >>>..
                  ''')

    def testSubjectExtendsLeft2(self):
        self.check('''
                                    ....>>>...........
                                      ..<<<...........
                             ...........>>>..
                  ''')

    def testSubjectExtendsRight1(self):
        self.check('''
                                        >>>....
                                        <<<..
                                      ..>>>
                  ''')

    def testSubjectExtendsBoth(self):
        self.check('''
                                  ......>>>...........
                                      ..<<<...
                                     ...>>>..
                  ''')

    def testQueryExtendsLeft1(self):
        self.check('''
                                    ....>>>
                                  ......<<<
                                        >>>......
                  ''')

    def testQueryExtendsLeft2(self):
        self.check('''
                                    ....>>>...........
                                  ......<<<...........
                             ...........>>>......
                  ''')

    def testQueryExtendsRight1(self):
        self.check('''
                                        >>>
                                        <<<..
                                      ..>>>
                  ''')

    def testQueryExtendsRight2(self):
        self.check('''
                                    ....>>>
                                    ....<<<..
                                      ..>>>....
                  ''')

    def testQueryExtendsBoth(self):
        self.check('''
                                    ....>>>...........
                                  ......<<<...............
                         ...............>>>......
                  ''')

    def testSubjectExtendsLeftQueryExtendsRight(self):
        self.check('''
                                    ....>>>...........
                                      ..<<<...............
                         ...............>>>..
                  ''')

    def testSubjectExtendsRightQueryExtendsLeft(self):
        self.check('''
                                      ..<<<...............
                                    ....>>>...........
                             ...........>>>....
                  ''')


class SubjectNegativeQueryPositive(TestCase, TestNormalizeHSPMixin):

    def testIdentical1(self):
        self.check('''
                                    <
                                    >
                                    >
                  ''')

    def testIdentical2(self):
        self.check('''
                                    <<<<
                                    >>>>
                                    >>>>
                  ''')

    def testSubjectExtendsLeft1(self):
        self.check('''
                                             ....>>>
                                               ..<<<
                                                 >>>..
                   ''')

    def testSubjectExtendsLeft2(self):
        self.check('''
                                             ....>>>...........
                                               ..<<<...........
                                      ...........>>>..
                   ''')

    def testSubjectExtendsRight1(self):
        self.check('''
                                             >>>.................
                                             <<<...........
                                  ...........>>>
                   ''')

    def testSubjectExtendsRight2(self):
        self.check('''
                                             ....>>>.................
                                             ....<<<...........
                                      ...........>>>....
                   ''')

    def testSubjectExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..>>>...
                                     ...>>>..
                  ''')

    def testQueryExtendsLeft1(self):
        self.check('''
                                               ..<<<
                                             ....>>>
                                                 >>>....
                   ''')

    def testQueryExtendsLeft2(self):
        self.check('''
                                               ..<<<...........
                                             ....>>>...........
                                      ...........>>>....
                   ''')

    def testQueryExtendsRight1(self):
        self.check('''
                                             <<<...........
                                             >>>.................
                            .................>>>
                   ''')

    def testQueryExtendsRight2(self):
        self.check('''
                                             ....<<<...........
                                             ....>>>.................
                                .................>>>....
                   ''')

    def testQueryExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..>>>...
                                     ...>>>..
                  ''')

    def testSubjectExtendsLeftQueryExtendsRight(self):
        self.check('''
                                    ....<<<...........
                                      ..>>>...............
                         ...............>>>..
                  ''')

    def testSubjectExtendsRightQueryExtendsLeft(self):
        self.check('''
                                      ..<<<...............
                                    ....>>>...........
                             ...........>>>....
                  ''')


class SubjectNegativeQueryNegative(TestCase, TestNormalizeHSPMixin):

    def testIdentical1(self):
        self.check('''
                                    <
                                    <
                                    >
                  ''')

    def testIdentical2(self):
        self.check('''
                                    <<<<
                                    <<<<
                                    >>>>
                  ''')

    def testSubjectExtendsLeft1(self):
        self.check('''
                                             ....<<<
                                               ..<<<
                                               ..>>>
                   ''')

    def testSubjectExtendsLeft2(self):
        self.check('''
                                             ....<<<...........
                                               ..<<<...........
                                               ..>>>...........
                   ''')

    def testSubjectExtendsRight1(self):
        self.check('''
                                             <<<.................
                                             <<<...........
                                             <<<...........
                   ''')

    def testSubjectExtendsRight2(self):
        self.check('''
                                             ....<<<.................
                                             ....<<<...........
                                             ....<<<...........
                   ''')

    def testSubjectExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..<<<...
                                      ..<<<...
                  ''')

    def testQueryExtendsLeft1(self):
        self.check('''
                                               ..<<<
                                             ....<<<
                                             ....<<<
                   ''')

    def testQueryExtendsLeft2(self):
        self.check('''
                                               ..<<<...........
                                             ....<<<...........
                                             ....<<<...........
                   ''')

    def testQueryExtendsRight1(self):
        self.check('''
                                             <<<...........
                                             <<<.................
                                             <<<.................
                   ''')

    def testQueryExtendsRight2(self):
        self.check('''
                                             ....<<<...........
                                             ....<<<.................
                                             ....<<<.................
                   ''')

    def testQueryExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..<<<...
                                      ..<<<...
                  ''')

    def testSubjectExtendsLeftQueryExtendsRight(self):
        self.check('''
                                    ....<<<...........
                                      ..<<<...............
                                      ..<<<...............
                  ''')

    def testSubjectExtendsRightQueryExtendsLeft(self):
        self.check('''
                                      ..<<<...............
                                    ....<<<...........
                                    ....<<<...........
                  ''')


class Old_QueryPositiveSubjectPositive(TestCase):
    """
    Tests for normalizeHSP when the subject start is less than the subject end.

    NOTE: Please don't add any tests below. Use the Template based tests above.
    """

    frame = Frame(query=1, subject=1)

    def testIdentical(self):
        """The subject start and end are identical to those of the query.

             ssss
             qqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=4,
                  frame=self.frame)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': 0,
            'queryEnd': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """The subject overlaps the query to the left.

              ssssss
                qqqq
        """
        hsp = HSP(subjectStart=3, subjectEnd=6, queryStart=1, queryEnd=4,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=6,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=4,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=5,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=2, subjectEnd=5, queryStart=1, queryEnd=4,
                  frame=self.frame)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'queryStart': 1,
            'queryEnd': 5,
        }, normalized)


class Old_QueryNegativeSubjectPositive(TestCase):
    """
    Tests for normalizeHSP when the query end is less than the query start.

    NOTE: Please don't add any tests below. Use the Template based tests above.
    """

    frame = Frame(query=-1, subject=1)

    def testIdentical(self):
        """The subject start and end are identical to those of the query.

             ssss
             qqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=4, queryEnd=1,
                  frame=self.frame)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': 0,
            'queryEnd': 4,
        }, normalized)

    def testSubjectExtendsRight(self):
        """The subject overlaps the query to the right.

              ssssss
              qqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=6, queryStart=4, queryEnd=1,
                  frame=self.frame)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 6,
            'queryStart': 0,
            'queryEnd': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """The subject overlaps the query to the left.

              ssssss
                qqqq
        """
        hsp = HSP(subjectStart=3, subjectEnd=6, queryStart=4, queryEnd=1,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=6, queryEnd=3,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=6, queryEnd=1,
                  frame=self.frame)
        normalized = normalizeHSP(hsp, 6)
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'queryStart': -2,
            'queryEnd': 4,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """The query extends to the right and left of the subject.

                ssss
               qqqqqq
        """
        hsp = HSP(subjectStart=1, subjectEnd=4, queryStart=5, queryEnd=2,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=2, subjectEnd=5, queryStart=4, queryEnd=1,
                  frame=self.frame)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'queryStart': 1,
            'queryEnd': 5,
        }, normalized)


class Old_QueryPositiveSubjectNegative(TestCase):
    """
    Tests for normalizeHSP when the subject start is greater than the subject
    end.

    NOTE: Please don't add any tests below. Use the Template based tests above.
    """
    frame = Frame(query=1, subject=-1)

    def testIdentical(self):
        """
        The subject start and end are identical to those of the query.

              ssss
              qqqq
        """
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=1, queryEnd=4,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=1, queryEnd=4,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=3, queryEnd=6,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=3, queryEnd=6,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=2, subjectEnd=1, queryStart=5, queryEnd=6,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=4, subjectEnd=1, queryStart=3, queryEnd=6,
                  frame=self.frame)
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
        hsp = HSP(subjectStart=5, subjectEnd=2, queryStart=1, queryEnd=4,
                  frame=self.frame)
        normalized = normalizeHSP(hsp, 4)
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'queryStart': 1,
            'queryEnd': 5,
        }, normalized)

    def test20130721Debugging(self):
        """
        This is an example I manually examined on 2013-07-21.
        """
        hsp = HSP(subjectStart=9018, subjectEnd=8764, queryStart=66,
                  queryEnd=315, frame=self.frame)
        normalized = normalizeHSP(hsp, 316)
        self.assertEqual({
            'subjectStart': 8763,
            'subjectEnd': 9018,
            'queryStart': 8767,
            'queryEnd': 9083,
        }, normalized)

    def test20131113Debugging(self):
        """
        This is an example I manually examined on 2013-11-13.
        """
        subject = (
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGTGTCCGGCA'
            'AGGTTGCCAAGGAGCAGATCGACATCGATAACGCCAAGCACACCAAGTGATGCACTGA'
            'CGACGGGTGAGGCCCAGATTCCTACGGCCTGGGCCTCTGTCTGCGTCGGGATGCCATT'
            'AGGCCGGTAGGATCGGTCACATGATCGATCCCAAGCTCCTGCGAACGGATCCGGACGC'
            'CGTTCGTCGCTCCCAGGCCGCCCGCGGCGAGGACTCCTCGGTTGTGGACGACGTTGTC'
            'GCCGCAGATGAGGCTCGTCGTGAGGCTATTGCTGCCCATGAGAACCTGCGTGCAGAAC'
            'AGAAGGGACTCGGCAAGCGAATCGCTAAAGCATCCGGTG')

        query = (
            'GTC-AGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGTGTCCGGCA'
            'AGGTTGCCAAGGAGCAGATCGACATCGATAACGCCAAGCACACCAAGTGATGCACTGA'
            'CGACGGGTGAGGCCCAGATTCCTACGGCCTGGGCCTCTGTCTGCGTCGGGATGCCATT'
            'AGGCCGCTAGGATCGGTCACATGATCGATCCCAAGCTCCTGCGAACGGATCCGGACGC'
            'CGTTCGTCGCTCCCAGGCCGCCCGCGGCGAGGACTCCTCGGTTGTGGACGACGTTGTC'
            'GCCGCAGATGAGGCTCGTCGTGAGGCTATTGCTGCCCATGAGAACCTGCGTGCAGAAC'
            'AGAAGGGACTCGGCAAGCGAATCGCTAAAGCATCCGGTG')

        hsp = HSP(subjectStart=2339751, subjectEnd=2339365, queryStart=1,
                  queryEnd=386, frame=self.frame, subject=subject, query=query)
        normalized = normalizeHSP(hsp, 396)
        self.assertEqual({
            'subjectStart': 2339364,
            'subjectEnd': 2339751,
            'queryStart': 2339355,
            'queryEnd': 2339751,
        }, normalized)

from unittest import TestCase

from dark.diamond.hsp import normalizeHSP


class FakeHSP(dict):
    def __init__(self, subjectStart, subjectEnd, queryStart, queryEnd, frame,
                 subject='', query='', bitscore=None, evalue=None, btop=''):
        """
        A fake HSP class with 1-based offsets and key names as used by DIAMOND.
        """
        if frame > 0:
            if not queryStart < queryEnd:
                raise ValueError('queryStart (%d) not less than queryEnd '
                                 '(%d) when frame (%d) is positive.' %
                                 (queryStart, queryEnd, frame))
        else:
            if not queryStart > queryEnd:
                raise ValueError('queryStart (%d) not greater than queryEnd '
                                 '(%d) when frame (%d) is negative.' %
                                 (queryStart, queryEnd, frame))

        self['sbjct_start'] = subjectStart
        self['sbjct_end'] = subjectEnd
        self['query_start'] = queryStart
        self['query_end'] = queryEnd
        self['frame'] = frame
        self['sbjct'] = subject
        self['query'] = query
        self['bits'] = bitscore
        self['expect'] = evalue
        self['btop'] = btop


class TestBlastxFramePlus1(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=1 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the first nucleotide.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects.

    def testIdentical(self):
        """
        The subject start and end are identical to those of the translated
        query.

             ....
             ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=12,
                      frame=1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=1, queryEnd=12,
                      frame=1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 2,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=7, queryEnd=18,
                      frame=1)
        normalized = normalizeHSP(hsp, 18, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=12,
                      frame=1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=12,
                      frame=1)
        normalized = normalizeHSP(hsp, 18, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=7, queryEnd=18,
                      frame=1)
        normalized = normalizeHSP(hsp, 21, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 5,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=5, queryStart=3, queryEnd=12,
                      frame=1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 1,
            'readEndInSubject': 5,
        }, normalized)


class TestBlastxFramePlus2(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=2 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the first nucleotide.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects.

    def testIdentical(self):
        """
        The subject start and end are identical to those of the translated
        query.

             ....
             ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=13,
                      frame=2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=2, queryEnd=13,
                      frame=2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 2,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=8, queryEnd=19,
                      frame=2)
        normalized = normalizeHSP(hsp, 19, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=13,
                      frame=2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=13,
                      frame=2)
        normalized = normalizeHSP(hsp, 19, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=8, queryEnd=19,
                      frame=2)
        normalized = normalizeHSP(hsp, 22, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 5,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=5, queryStart=4, queryEnd=13,
                      frame=2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 1,
            'readEndInSubject': 5,
        }, normalized)


class TestBlastxFramePlus3(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=3 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the first nucleotide.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects.

    def testIdentical(self):
        """
        The subject start and end are identical to those of the translated
        query.

             ....
             ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=14,
                      frame=3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=3, queryEnd=14,
                      frame=3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 2,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=9, queryEnd=20,
                      frame=3)
        normalized = normalizeHSP(hsp, 20, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=14,
                      frame=3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=14,
                      frame=3)
        normalized = normalizeHSP(hsp, 20, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=9, queryEnd=20,
                      frame=3)
        normalized = normalizeHSP(hsp, 23, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 5,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=5, queryStart=5, queryEnd=14,
                      frame=3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 1,
            'readEndInSubject': 5,
        }, normalized)


class TestBlastxFrameMinus1(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=-1 (i.e., the
    query matches in the reverse (complemented) order it was given to DIAMOND,
    and the translation frame starts at the last nucleotide.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects.

    def testIdentical(self):
        """
        If the query and subject match completely, the normalized HSP must
        be correct.

               ....
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=12, queryEnd=1,
                      frame=-1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 2,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-1)
        normalized = normalizeHSP(hsp, 18, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=18, queryEnd=7,
                      frame=-1)
        normalized = normalizeHSP(hsp, 18, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=15, queryEnd=4,
                      frame=-1)
        normalized = normalizeHSP(hsp, 21, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 5,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=4, queryStart=12, queryEnd=4,
                      frame=-1)
        normalized = normalizeHSP(hsp, 12, 'blastx')
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 3,
            'readStartInSubject': 1,
            'readEndInSubject': 5,
        }, normalized)


class TestBlastxFrameMinus2(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=-2 (i.e., the
    query matches in the reverse (complemented) order it was given to DIAMOND,
    and the translation frame starts at the last nucleotide.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects.

    def testIdentical(self):
        """
        If the query and subject match completely, the normalized HSP must
        be correct.

               ....
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=12, queryEnd=1,
                      frame=-2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 2,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-2)
        normalized = normalizeHSP(hsp, 19, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=18, queryEnd=7,
                      frame=-2)
        normalized = normalizeHSP(hsp, 19, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=15, queryEnd=4,
                      frame=-2)
        normalized = normalizeHSP(hsp, 22, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 5,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=4, queryStart=12, queryEnd=4,
                      frame=-2)
        normalized = normalizeHSP(hsp, 13, 'blastx')
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 3,
            'readStartInSubject': 1,
            'readEndInSubject': 5,
        }, normalized)


class TestBlastxFrameMinus3(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=-3 (i.e., the
    query matches in the reverse (complemented) order it was given to DIAMOND,
    and the translation frame starts at the last nucleotide.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects.

    def testIdentical(self):
        """
        If the query and subject match completely, the normalized HSP must
        be correct.

               ....
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=12, queryEnd=1,
                      frame=-3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 2,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-3)
        normalized = normalizeHSP(hsp, 20, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 4,
        }, normalized)

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1,
                      frame=-3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=18, queryEnd=7,
                      frame=-3)
        normalized = normalizeHSP(hsp, 20, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInSubject': 0,
            'readEndInSubject': 6,
        }, normalized)

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=15, queryEnd=4,
                      frame=-3)
        normalized = normalizeHSP(hsp, 23, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInSubject': -2,
            'readEndInSubject': 5,
        }, normalized)

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=4, queryStart=12, queryEnd=4,
                      frame=-3)
        normalized = normalizeHSP(hsp, 14, 'blastx')
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 4,
            'readStart': 0,
            'readEnd': 3,
            'readStartInSubject': 1,
            'readEndInSubject': 5,
        }, normalized)

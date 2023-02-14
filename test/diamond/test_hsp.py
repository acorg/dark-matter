from unittest import TestCase

from dark.diamond.hsp import normalizeHSP


class FakeHSP(dict):
    def __init__(
        self,
        subjectStart,
        subjectEnd,
        queryStart,
        queryEnd,
        frame,
        subject="",
        query="",
        bitscore=None,
        evalue=None,
        btop="",
    ):
        """
        A fake HSP class with 1-based offsets and key names as used by DIAMOND.
        """
        if frame > 0:
            if not queryStart < queryEnd:
                raise ValueError(
                    "queryStart (%d) not less than queryEnd "
                    "(%d) when frame (%d) is positive." % (queryStart, queryEnd, frame)
                )
        else:
            if not queryStart > queryEnd:
                raise ValueError(
                    "queryStart (%d) not greater than queryEnd "
                    "(%d) when frame (%d) is negative." % (queryStart, queryEnd, frame)
                )

        self["sbjct_start"] = subjectStart
        self["sbjct_end"] = subjectEnd
        self["query_start"] = queryStart
        self["query_end"] = queryEnd
        self["frame"] = frame
        self["sbjct"] = subject
        self["query"] = query
        self["bits"] = bitscore
        self["expect"] = evalue
        self["btop"] = btop


class TestBlastxFramePlus1(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=1 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the first nucleotide).
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
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=12, frame=1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=1, queryEnd=12, frame=1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=7, queryEnd=18, frame=1)
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=12, frame=1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=12, frame=1)
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=7, queryEnd=18, frame=1)
        normalized = normalizeHSP(hsp, 21, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=5, queryStart=3, queryEnd=12, frame=1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFramePlus2(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=2 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the second nucleotide).
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
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=13, frame=2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=2, queryEnd=13, frame=2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=8, queryEnd=19, frame=2)
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=13, frame=2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=2, queryEnd=13, frame=2)
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=8, queryEnd=19, frame=2)
        normalized = normalizeHSP(hsp, 22, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=5, queryStart=4, queryEnd=13, frame=2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFramePlus3(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=3 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the third nucleotide).
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
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=14, frame=3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=3, queryEnd=14, frame=3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=9, queryEnd=20, frame=3)
        normalized = normalizeHSP(hsp, 20, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=14, frame=3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=3, queryEnd=14, frame=3)
        normalized = normalizeHSP(hsp, 20, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=9, queryEnd=20, frame=3)
        normalized = normalizeHSP(hsp, 23, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=5, queryStart=5, queryEnd=14, frame=3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFrameMinus1(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=-1 (i.e., the
    query matches in the reverse (complemented) order it was given to DIAMOND,
    and the translation frame starts at the last nucleotide).
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
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=12, queryEnd=1, frame=-1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-1)
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=18, queryEnd=7, frame=-1)
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=15, queryEnd=4, frame=-1)
        normalized = normalizeHSP(hsp, 21, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=4, queryStart=12, queryEnd=4, frame=-1)
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFrameMinus2(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=-2 (i.e., the
    query matches in the reverse (complemented) order it was given to DIAMOND,
    and the translation frame starts at the second last nucleotide).
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
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=12, queryEnd=1, frame=-2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-2)
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=18, queryEnd=7, frame=-2)
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=15, queryEnd=4, frame=-2)
        normalized = normalizeHSP(hsp, 22, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=4, queryStart=12, queryEnd=4, frame=-2)
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFrameMinus3(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=-3 (i.e., the
    query matches in the reverse (complemented) order it was given to DIAMOND,
    and the translation frame starts at the third last nucleotide).
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
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeft(self):
        """
        The subject overlaps the translated query to the left.

              ss....
                ....
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, queryStart=12, queryEnd=1, frame=-3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeft(self):
        """
        The translated query extends to the left of the subject.

               ....
             qq....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-3)
        normalized = normalizeHSP(hsp, 20, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRight(self):
        """
        The subject extends to the right of the translated query.

               ....ss
               ....
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=12, queryEnd=1, frame=-3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRight(self):
        """
        The translated query extends to the right of the subject.

               ....
               ....qq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=18, queryEnd=7, frame=-3)
        normalized = normalizeHSP(hsp, 20, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeft(self):
        """
        The translated query extends to the right and left of the subject.

                ....
              qq....q
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, queryStart=15, queryEnd=4, frame=-3)
        normalized = normalizeHSP(hsp, 23, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeft(self):
        """
        The subject extends to the right and left of the translated query.

               s...sss
                ...q
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=4, queryStart=12, queryEnd=4, frame=-3)
        normalized = normalizeHSP(hsp, 14, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFramePlus1WithGaps(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=1 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the first nucleotide) and with gaps.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects. Gaps are shown as a hyphen.

    def testIdenticalWithQueryGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the query.

             ....
             .-..
        """
        hsp = FakeHSP(
            subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=9, frame=1, btop="1-K2"
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithTwoQueryGaps(self):
        """
        The subject start and end are identical to those of the translated
        query, given two gaps in the query.

             ....
             .--.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=1,
            queryEnd=6,
            frame=1,
            btop="1-K-K1",
        )
        normalized = normalizeHSP(hsp, 6, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 2,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithSubjectGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the subject.

             .-..
             ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=1,
            queryEnd=12,
            frame=1,
            btop="1K-2",
        )
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithTwoSubjectGaps(self):
        """
        The subject start and end are identical to those of the translated
        query, given two gaps in the subject.

             .--.
             ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=2,
            queryStart=1,
            queryEnd=12,
            frame=1,
            btop="1K-K-1",
        )
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 2,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithQueryAndSubjectGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the query and one in the subject.

             .-..
             ..-.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=1,
            queryEnd=9,
            frame=1,
            btop="1K--K2",
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeftWithQueryGap(self):
        """
        The subject overlaps the translated query to the left. The query has a
        gap.

              ss....
                ..-.
        """
        hsp = FakeHSP(
            subjectStart=3, subjectEnd=6, queryStart=1, queryEnd=9, frame=1, btop="2-K1"
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testSubjectExtendsLeftWithSubjectGap(self):
        """
        The subject overlaps the translated query to the left. The subject has
        a gap.

              ss.-..
                ....
        """
        hsp = FakeHSP(
            subjectStart=3,
            subjectEnd=5,
            queryStart=1,
            queryEnd=12,
            frame=1,
            btop="1K-2",
        )
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testSubjectExtendsLeftWithQueryAndSubjectGap(self):
        """
        The subject overlaps the translated query to the left. The query and
        subject both have a gap.

              ss.-..
                ..-.
        """
        hsp = FakeHSP(
            subjectStart=3,
            subjectEnd=5,
            queryStart=1,
            queryEnd=9,
            frame=1,
            btop="1K--K1",
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeftWithQueryGap(self):
        """
        The translated query extends to the left of the subject and the query
        has a gap.

               ....
             qq..-.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=7,
            queryEnd=15,
            frame=1,
            btop="2-K1",
        )
        normalized = normalizeHSP(hsp, 15, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsLeftWithSubjectGap(self):
        """
        The translated query extends to the left of the subject and the subject
        has a gap.

               ..-.
             qq....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=7,
            queryEnd=18,
            frame=1,
            btop="2K-1",
        )
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsLeftWithQueryAndSubjectGap(self):
        """
        The translated query extends to the left of the subject and the query
        and subject have gaps.

               ..-.
             qq.-..
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=7,
            queryEnd=15,
            frame=1,
            btop="1-KK-1",
        )
        normalized = normalizeHSP(hsp, 15, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRightWithQueryGap(self):
        """
        The subject extends to the right of the translated query and the query
        has a gap.

               ....ss
               ..-.
        """
        hsp = FakeHSP(
            subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=9, frame=1, btop="2-K1"
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRightWithSubjectGap(self):
        """
        The subject extends to the right of the translated query and the
        subject has a gap.

               ..-.ss
               ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=1,
            queryEnd=12,
            frame=1,
            btop="2K-1",
        )
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRightWithQueryAndSubjectGap(self):
        """
        The subject extends to the right of the translated query and the
        query and subject both have gaps.

               ..-.ss
               .-..
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=1,
            queryEnd=9,
            frame=1,
            btop="2-KK-1",
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRightWithQueryGap(self):
        """
        The translated query extends to the right of the subject and the
        query has a gap.

               ....
               ..-.qq
        """
        hsp = FakeHSP(
            subjectStart=1, subjectEnd=4, queryStart=1, queryEnd=9, frame=1, btop="2-K1"
        )
        normalized = normalizeHSP(hsp, 15, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightWithSubjectGap(self):
        """
        The translated query extends to the right of the subject and the
        subject has a gap.

               ..-.
               ....qq
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=1,
            queryEnd=12,
            frame=1,
            btop="2K-1",
        )
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightWithQueryAndSubjectGap(self):
        """
        The translated query extends to the right of the subject and the
        query and subject both have gaps.

               ..-.
               .-..qq
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=1,
            queryEnd=9,
            frame=1,
            btop="1-KK-1",
        )
        normalized = normalizeHSP(hsp, 15, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeftWithQueryGap(self):
        """
        The translated query extends to the right and left of the subject and
        the query has a gap.

                ....
              qq.-..q
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=7,
            queryEnd=15,
            frame=1,
            btop="1-K2",
        )
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeftWithSubjectGap(self):
        """
        The translated query extends to the right and left of the subject and
        the subject has a gap.

                .-..
              qq....q
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=7,
            queryEnd=18,
            frame=1,
            btop="1K-2",
        )
        normalized = normalizeHSP(hsp, 21, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeftWithQueryAndSubjectGap(self):
        """
        The translated query extends to the right and left of the subject and
        the query and the subject have gaps.

                .-..
              qq..-.q
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=7,
            queryEnd=15,
            frame=1,
            btop="1K--K1",
        )
        normalized = normalizeHSP(hsp, 18, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeftWithQueryGap(self):
        """
        The subject extends to the right and left of the translated query and
        the query has a gap.

               s....ss
                .-..
        """
        hsp = FakeHSP(
            subjectStart=2, subjectEnd=5, queryStart=3, queryEnd=9, frame=1, btop="1-K1"
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeftWithSubjectGap(self):
        """
        The subject extends to the right and left of the translated query and
        the subject has a gap.

               s.-..ss
                ....
        """
        hsp = FakeHSP(
            subjectStart=2,
            subjectEnd=4,
            queryStart=3,
            queryEnd=12,
            frame=1,
            btop="1K-1",
        )
        normalized = normalizeHSP(hsp, 12, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeftWithQueryAndSubjectGap(self):
        """
        The subject extends to the right and left of the translated query and
        the query and subjects have gaps.

               s.-..ss
                ..-.
        """
        hsp = FakeHSP(
            subjectStart=2,
            subjectEnd=4,
            queryStart=3,
            queryEnd=9,
            frame=1,
            btop="1K--K1",
        )
        normalized = normalizeHSP(hsp, 9, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFramePlus2WithGaps(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=2 (i.e., the
    query matches in the order it was given to DIAMOND, and the translation
    frame starts at the second nucleotide) and with gaps.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects. Gaps are shown as a hyphen.

    def testIdenticalWithQueryGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the query.

             ....
             .-..
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="1-K2",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithTwoQueryGaps(self):
        """
        The subject start and end are identical to those of the translated
        query, given two gaps in the query.

             ....
             .--.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=2,
            queryEnd=7,
            frame=2,
            btop="1-K-K1",
        )
        normalized = normalizeHSP(hsp, 7, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 2,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithSubjectGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the subject.

             .-..
             ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=2,
            queryEnd=13,
            frame=2,
            btop="1K-2",
        )
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithTwoSubjectGaps(self):
        """
        The subject start and end are identical to those of the translated
        query, given two gaps in the subject.

             .--.
             ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=2,
            queryStart=2,
            queryEnd=13,
            frame=2,
            btop="1K-K-1",
        )
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 2,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithQueryAndSubjectGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the query and one in the subject.

             .-..
             ..-.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="1K--K2",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsLeftWithQueryGap(self):
        """
        The subject overlaps the translated query to the left. The query has a
        gap.

              ss....
                ..-.
        """
        hsp = FakeHSP(
            subjectStart=3,
            subjectEnd=6,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="2-K1",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 6,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testSubjectExtendsLeftWithSubjectGap(self):
        """
        The subject overlaps the translated query to the left. The subject has
        a gap.

              ss.-..
                ....
        """
        hsp = FakeHSP(
            subjectStart=3,
            subjectEnd=5,
            queryStart=2,
            queryEnd=13,
            frame=2,
            btop="1K-2",
        )
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testSubjectExtendsLeftWithQueryAndSubjectGap(self):
        """
        The subject overlaps the translated query to the left. The query and
        subject both have a gap.

              ss.-..
                ..-.
        """
        hsp = FakeHSP(
            subjectStart=3,
            subjectEnd=5,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="1K--K1",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 2,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 2,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsLeftWithQueryGap(self):
        """
        The translated query extends to the left of the subject and the query
        has a gap.

               ....
             qq..-.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=8,
            queryEnd=16,
            frame=2,
            btop="2-K1",
        )
        normalized = normalizeHSP(hsp, 16, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsLeftWithSubjectGap(self):
        """
        The translated query extends to the left of the subject and the subject
        has a gap.

               ..-.
             qq....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=8,
            queryEnd=19,
            frame=2,
            btop="2K-1",
        )
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsLeftWithQueryAndSubjectGap(self):
        """
        The translated query extends to the left of the subject and the query
        and subject have gaps.

               ..-.
             qq.-..
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=8,
            queryEnd=16,
            frame=2,
            btop="1-KK-1",
        )
        normalized = normalizeHSP(hsp, 16, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRightWithQueryGap(self):
        """
        The subject extends to the right of the translated query and the query
        has a gap.

               ....ss
               ..-.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="2-K1",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRightWithSubjectGap(self):
        """
        The subject extends to the right of the translated query and the
        subject has a gap.

               ..-.ss
               ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=2,
            queryEnd=13,
            frame=2,
            btop="2K-1",
        )
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testSubjectExtendsRightWithQueryAndSubjectGap(self):
        """
        The subject extends to the right of the translated query and the
        query and subject both have gaps.

               ..-.ss
               .-..
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="2-KK-1",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testQueryExtendsRightWithQueryGap(self):
        """
        The translated query extends to the right of the subject and the
        query has a gap.

               ....
               ..-.qq
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="2-K1",
        )
        normalized = normalizeHSP(hsp, 16, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightWithSubjectGap(self):
        """
        The translated query extends to the right of the subject and the
        subject has a gap.

               ..-.
               ....qq
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=2,
            queryEnd=13,
            frame=2,
            btop="2K-1",
        )
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightWithQueryAndSubjectGap(self):
        """
        The translated query extends to the right of the subject and the
        query and subject both have gaps.

               ..-.
               .-..qq
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=2,
            queryEnd=10,
            frame=2,
            btop="1-KK-1",
        )
        normalized = normalizeHSP(hsp, 16, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 6,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeftWithQueryGap(self):
        """
        The translated query extends to the right and left of the subject and
        the query has a gap.

                ....
              qq.-..q
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=8,
            queryEnd=16,
            frame=2,
            btop="1-K2",
        )
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeftWithSubjectGap(self):
        """
        The translated query extends to the right and left of the subject and
        the subject has a gap.

                .-..
              qq....q
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=8,
            queryEnd=19,
            frame=2,
            btop="1K-2",
        )
        normalized = normalizeHSP(hsp, 22, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 6,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testQueryExtendsRightAndLeftWithQueryAndSubjectGap(self):
        """
        The translated query extends to the right and left of the subject and
        the query and the subject have gaps.

                .-..
              qq..-.q
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=8,
            queryEnd=16,
            frame=2,
            btop="1K--K1",
        )
        normalized = normalizeHSP(hsp, 19, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 2,
                "readEnd": 5,
                "readStartInSubject": -2,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeftWithQueryGap(self):
        """
        The subject extends to the right and left of the translated query and
        the query has a gap.

               s....ss
                .-..
        """
        hsp = FakeHSP(
            subjectStart=2,
            subjectEnd=5,
            queryStart=4,
            queryEnd=10,
            frame=2,
            btop="1-K1",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 5,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeftWithSubjectGap(self):
        """
        The subject extends to the right and left of the translated query and
        the subject has a gap.

               s.-..ss
                ....
        """
        hsp = FakeHSP(
            subjectStart=2,
            subjectEnd=4,
            queryStart=4,
            queryEnd=13,
            frame=2,
            btop="1K-1",
        )
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )

    def testSubjectExtendsRightAndLeftWithQueryAndSubjectGap(self):
        """
        The subject extends to the right and left of the translated query and
        the query and subjects have gaps.

               s.-..ss
                ..-.
        """
        hsp = FakeHSP(
            subjectStart=2,
            subjectEnd=4,
            queryStart=4,
            queryEnd=10,
            frame=2,
            btop="1K--K1",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 1,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 1,
                "readEndInSubject": 5,
            },
            normalized,
        )


class TestBlastxFrameMinus1WithGaps(TestCase):
    """
    Tests for normalizeHSP for DIAMOND blastx output when frame=-1 (i.e., the
    query matches in the reverse (complemented) order it was given to DIAMOND,
    and the translation frame starts at the last nucleotide) and with gaps.
    """

    # All query offsets and lengths must be in terms of nucleotides.
    # Subject offsets are in terms of protein AA sequences.  This is how
    # DIAMOND reports those offsets.
    #
    # In the little diagrams in the docstrings below, the first line is the
    # subject and the second the query. Dots indicate where the matched
    # region is. The queries are shown translated so as to line up properly
    # with the subjects. Gaps are shown as a hyphen.

    def testIdenticalWithQueryGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the query.

             ....
             .-..
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=10,
            queryEnd=2,
            frame=-1,
            btop="1-K2",
        )
        normalized = normalizeHSP(hsp, 10, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 3,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithTwoQueryGaps(self):
        """
        The subject start and end are identical to those of the translated
        query, given two gaps in the query.

             ....
             .--.
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=4,
            queryStart=7,
            queryEnd=2,
            frame=-1,
            btop="1-K-K1",
        )
        normalized = normalizeHSP(hsp, 7, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 4,
                "readStart": 0,
                "readEnd": 2,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithSubjectGap(self):
        """
        The subject start and end are identical to those of the translated
        query, given one gap in the subject.

             .-..
             ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=3,
            queryStart=13,
            queryEnd=2,
            frame=-1,
            btop="1K-2",
        )
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 3,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

    def testIdenticalWithTwoSubjectGaps(self):
        """
        The subject start and end are identical to those of the translated
        query, given two gaps in the subject.

             .--.
             ....
        """
        hsp = FakeHSP(
            subjectStart=1,
            subjectEnd=2,
            queryStart=13,
            queryEnd=2,
            frame=-1,
            btop="1K-K-1",
        )
        normalized = normalizeHSP(hsp, 13, "blastx")
        self.assertEqual(
            {
                "subjectStart": 0,
                "subjectEnd": 2,
                "readStart": 0,
                "readEnd": 4,
                "readStartInSubject": 0,
                "readEndInSubject": 4,
            },
            normalized,
        )

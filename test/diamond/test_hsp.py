from unittest import TestCase

from dark.diamond.hsp import normalizeHSP


class Frame(object):
    def __init__(self, read):
        self.read = read


class FakeHSP(dict):
    def __init__(self, subjectStart, subjectEnd, readStart, readEnd, frame,
                 hit='', read=''):
        """
        A fake HSP class (with 1-based offsets, as are used in DIAMOND).
        """
        self['sbjct_start'] = subjectStart
        self['sbjct_end'] = subjectEnd
        self['query_start'] = readStart
        self['query_end'] = readEnd
        self['frame'] = frame.read
        self['sbjct'] = hit
        self['query'] = read

        # In case you're thinking of adding it, the following assertion is
        # not valid:
        #
        #   assert abs(subjectEnd - subjectStart) == abs(readEnd - readStart)
        #
        # That's because DIAMOND might find a match that requires a gap in
        # the read or in the hit. The indices that it reports do not
        # include the gap and so the differences in the lengths of the
        # sections of the read and hit may not be the same.


class Old_ReadPositiveHitPositive(TestCase):
    """
    Tests for normalizeHSP when the hit start is less than the hit end.

    NOTE: Please don't add any tests below. Use the Template based tests above.
    """

    frame = Frame(read=1)

    def testIdentical(self):
        """
        The hit start and end are identical to those of the read.

             ssss
             qqqq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, readStart=3, readEnd=12,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 12,
            'readStartInSubject': 0,
            'readEndInSubject': 4,
        }, normalized)

    def testHitExtendsLeft(self):
        """
        The hit overlaps the read to the left.

              ssssss
                qqqq
        """
        hsp = FakeHSP(subjectStart=3, subjectEnd=6, readStart=3, readEnd=12,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastx')
        self.assertEqual({
            'subjectStart': 2,
            'subjectEnd': 6,
            'readStart': 2,
            'readEnd': 12,
            'readStartInSubject': 2,
            'readEndInSubject': 6,
        }, normalized)

    def testReadExtendsLeft(self):
        """
        The read sticks out to the left of the hit.

               ssss
             qqqqqq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, readStart=9, readEnd=18,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 8,
            'readEnd': 18,
            'readStartInSubject': -2,
            'readEndInSubject': 4,
        }, normalized)

    def testReadExtendsRight(self):
        """
        The read sticks out to the right of the hit.

               ssss
               qqqqqq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, readStart=3, readEnd=12,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 2,
            'readEnd': 12,
            'readStartInSubject': 0,
            'readEndInSubject': 6,
        }, normalized)

    def testReadExtendsRightAndLeft(self):
        """
        The read extends to the right and left of the hit.

                ssss
               qqqqqq
        """
        hsp = FakeHSP(subjectStart=1, subjectEnd=4, readStart=6, readEnd=15,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastx')
        self.assertEqual({
            'subjectStart': 0,
            'subjectEnd': 4,
            'readStart': 5,
            'readEnd': 15,
            'readStartInSubject': -1,
            'readEndInSubject': 5,
        }, normalized)

    def testHitExtendsRightAndLeft(self):
        """
        The hit extends to the right and left of the read.

               sssssss
                qqqq
        """
        hsp = FakeHSP(subjectStart=2, subjectEnd=5, readStart=3, readEnd=12,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastx')
        self.assertEqual({
            'subjectStart': 1,
            'subjectEnd': 5,
            'readStart': 2,
            'readEnd': 12,
            'readStartInSubject': 1,
            'readEndInSubject': 5,
        }, normalized)

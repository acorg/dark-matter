from unittest import TestCase

from dark.blast.hsp import normalizeHSP, EValueHSP


class TestEValueHSP(TestCase):
    """
    Tests of the L{dark.hsp.blast.EValueHSP} class.
    """

    def testExpectedAttributes(self):
        """
        An EValueHSP must have the expected attributes.
        """
        hsp = EValueHSP(readStart=1, readEnd=2,
                        readStartInHit=3, readEndInHit=4,
                        hitStart=5, hitEnd=6,
                        readMatchedSequence='aaa', hitMatchedSequence='ccc',
                        score=7)
        self.assertEqual(1, hsp.readStart)
        self.assertEqual(2, hsp.readEnd)
        self.assertEqual(3, hsp.readStartInHit)
        self.assertEqual(4, hsp.readEndInHit)
        self.assertEqual(5, hsp.hitStart)
        self.assertEqual(6, hsp.hitEnd)
        self.assertEqual('aaa', hsp.readMatchedSequence)
        self.assertEqual('ccc', hsp.hitMatchedSequence)
        self.assertEqual(7, hsp.score)

    def testEq(self):
        """
        Two HSPs must compare properly with ==
        """
        hsp1 = EValueHSP(score=7)
        hsp2 = EValueHSP(score=7)
        self.assertEqual(hsp1, hsp2)

    def testLt(self):
        """
        Two HSPs must compare properly with <
        """
        hsp1 = EValueHSP(score=8)
        hsp2 = EValueHSP(score=7)
        self.assertTrue(hsp1 < hsp2)

    def testBetterThanTrue(self):
        """
        Two HSP scores must compare properly with betterThan if the passed
        score is worse than the score of the HSP.
        """
        self.assertTrue(EValueHSP(score=5).betterThan(7))

    def testBetterThanFalse(self):
        """
        Two HSP scores must compare properly with betterThan if the passed
        score is better than the score of the HSP.
        """
        self.assertFalse(EValueHSP(score=7).betterThan(5))


class Frame(object):
    def __init__(self, read, hit):
        self.read = read
        self.hit = hit


class FakeHSP(dict):
    def __init__(self, hitStart, hitEnd, readStart, readEnd, frame,
                 hit='', read=''):
        """
        A fake HSP class (with 1-based offsets, as are used in BLAST).
        """
        self['sbjct_start'] = hitStart
        self['sbjct_end'] = hitEnd
        self['query_start'] = readStart
        self['query_end'] = readEnd
        self['frame'] = (frame.read, frame.hit)
        self['sbjct'] = hit
        self['query'] = read

        # In case you're thinking of adding it, the following assertion is
        # not valid:
        #
        #   assert abs(hitEnd - hitStart) == abs(readEnd - readStart)
        #
        # That's because BLAST might find a match that requires a gap in
        # the read or in the hit. The indices that it reports do not
        # include the gap and so the differences in the lengths of the
        # sections of the read and hit may not be the same.


class Template(object):
    def __init__(self, template):
        template = template.split('\n')
        # Allow the first template line to be empty.
        if len(template[0]) == 0:
            template = template[1:]

        # Analyze the template hit.
        self.hit = template[0].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen, positive) = \
            self._analyze(self.hit)
        origin = spacesLen
        self.matchLen = matchLen
        self.hitLen = len(self.hit) - spacesLen
        self.hitMatchStart = leadingDotsLen
        self.hitPositive = positive

        # Analyze the template read.
        self.read = template[1].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen, positive) = \
            self._analyze(self.read)
        assert self.matchLen == matchLen
        self.readLen = len(self.read) - spacesLen
        self.readMatchStart = leadingDotsLen
        self.readPositive = positive

        # Analyze the template read result.
        self.readResult = template[2].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen, positive) = \
            self._analyze(self.readResult)
        assert self.matchLen == matchLen
        self.readResultStart = spacesLen - origin
        self.readResultLen = len(self.readResult) - spacesLen
        assert self.readResultLen == self.readLen
        self.readResultMatchStart = leadingDotsLen

    def leadingCharMatchLen(self, str, chars=' '):
        return len(str) - len(str.lstrip(chars))

    def _analyze(self, str):
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
        return FakeHSP(hitStart=self.hitMatchStart + 1,
                       hitEnd=self.hitMatchStart + self.matchLen,
                       readStart=self.readMatchStart + 1,
                       readEnd=self.readMatchStart + self.matchLen,
                       frame=Frame(
                           1 if self.readPositive else -1,
                           1 if self.hitPositive else -1,
                       ),
                       # Make non-random non-gapped read and hits.
                       hit='a' * self.matchLen,
                       read='a' * self.matchLen)


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

        self.assertTrue(template.hitPositive)
        self.assertEqual(18, template.hitLen)
        self.assertEqual(4, template.hitMatchStart)

        self.assertFalse(template.readPositive)
        self.assertEqual(20, template.readLen)
        self.assertEqual(2, template.readMatchStart)

        self.assertEqual(20, template.readResultLen)
        self.assertEqual(-11, template.readResultStart)
        self.assertEqual(15, template.readResultMatchStart)

        hsp = template.hsp()
        self.assertEqual(5, hsp['sbjct_start'])
        self.assertEqual(7, hsp['sbjct_end'])
        self.assertEqual(3, hsp['query_start'])
        self.assertEqual(5, hsp['query_end'])
        self.assertEqual((-1, 1), hsp['frame'])


class TestNormalizeHSPMixin(object):
    """
    Tests for normalizeHSP when the read and hit are both positive.
    """

    def check(self, templateStr):
        template = Template(templateStr)
        normalized = normalizeHSP(template.hsp(), template.readLen, 'blastn')
        self.assertEqual({
            'hitStart': template.hitMatchStart,
            'hitEnd': template.hitMatchStart + template.matchLen,
            'readStart': template.readMatchStart,
            'readEnd': template.readMatchStart + template.matchLen,
            'readStartInHit': template.readResultStart,
            'readEndInHit': template.readResultStart + template.readLen,
        }, normalized)


class HitPositiveReadPositive(TestCase, TestNormalizeHSPMixin):
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

    def testHitExtendsLeft1(self):
        self.check('''
                                    ....>>
                                        >>
                                        >>
                  ''')

    def testHitExtendsLeft2(self):
        self.check('''
                                    ....>>>...
                                        >>>...
                                        >>>...
                  ''')

    def testHitExtendsRight1(self):
        self.check('''
                                        >>......
                                        >>...
                                        >>...
                  ''')

    def testHitExtendsRight2(self):
        self.check('''
                                        ..>>>>......
                                        ..>>>>...
                                        ..>>>>...
                  ''')

    def testHitExtendsBoth(self):
        self.check('''
                                    ....>>>...........
                                      ..>>>....
                                      ..>>>....
                  ''')

    def testReadExtendsLeft1(self):
        self.check('''
                                        >>
                                    ....>>
                                    ....>>
                  ''')

    def testReadExtendsLeft2(self):
        self.check('''
                                        >>>...
                                    ....>>>...
                                    ....>>>...
                  ''')

    def testReadExtendsRight1(self):
        self.check('''
                                        >>...
                                        >>......
                                        >>......
                  ''')

    def testReadExtendsRight2(self):
        self.check('''
                                        ..>>>>...
                                        ..>>>>......
                                        ..>>>>......
                  ''')

    def testReadExtendsBoth(self):
        self.check('''
                                      ..>>>....
                                    ....>>>...........
                                    ....>>>...........
                  ''')

    def testHitExtendsLeftReadExtendsRight(self):
        self.check('''
                                    ....>>>...........
                                      ..>>>...............
                                      ..>>>...............
                  ''')

    def testHitExtendsRightReadExtendsLeft(self):
        self.check('''
                                      ..>>>...............
                                    ....>>>...........
                                    ....>>>...........
                  ''')


class HitPositiveReadNegative(TestCase, TestNormalizeHSPMixin):

    """
    This class appears to not be needed. As far as we have seen, the
    read is always positive with ascending start, end offsets.
    """

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

    def testHitExtendsLeft1(self):
        self.check('''
                                    ....>>>
                                      ..<<<
                                        >>>..
                  ''')

    def testHitExtendsLeft2(self):
        self.check('''
                                    ....>>>...........
                                      ..<<<...........
                             ...........>>>..
                  ''')

    def testHitExtendsRight1(self):
        self.check('''
                                        >>>....
                                        <<<..
                                      ..>>>
                  ''')

    def testHitExtendsBoth(self):
        self.check('''
                                  ......>>>...........
                                      ..<<<...
                                     ...>>>..
                  ''')

    def testReadExtendsLeft1(self):
        self.check('''
                                    ....>>>
                                  ......<<<
                                        >>>......
                  ''')

    def testReadExtendsLeft2(self):
        self.check('''
                                    ....>>>...........
                                  ......<<<...........
                             ...........>>>......
                  ''')

    def testReadExtendsRight1(self):
        self.check('''
                                        >>>
                                        <<<..
                                      ..>>>
                  ''')

    def testReadExtendsRight2(self):
        self.check('''
                                    ....>>>
                                    ....<<<..
                                      ..>>>....
                  ''')

    def testReadExtendsBoth(self):
        self.check('''
                                    ....>>>...........
                                  ......<<<...............
                         ...............>>>......
                  ''')

    def testHitExtendsLeftReadExtendsRight(self):
        self.check('''
                                    ....>>>...........
                                      ..<<<...............
                         ...............>>>..
                  ''')

    def testHitExtendsRightReadExtendsLeft(self):
        self.check('''
                                      ..<<<...............
                                    ....>>>...........
                             ...........>>>....
                  ''')


class HitNegativeReadPositive(TestCase, TestNormalizeHSPMixin):

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

    def testHitExtendsLeft1(self):
        self.check('''
                                             ....<<<
                                               ..>>>
                                                 >>>..
                   ''')

    def testHitExtendsLeft2(self):
        self.check('''
                                             ....<<<...........
                                               ..>>>...........
                                      ...........>>>..
                   ''')

    def testHitExtendsRight1(self):
        self.check('''
                                             <<<.................
                                             >>>...........
                                  ...........>>>
                   ''')

    def testHitExtendsRight2(self):
        self.check('''
                                             ....<<<.................
                                             ....>>>...........
                                      ...........>>>....
                   ''')

    def testHitExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..>>>...
                                     ...>>>..
                  ''')

    def testReadExtendsLeft1(self):
        self.check('''
                                               ..<<<
                                             ....>>>
                                                 >>>....
                   ''')

    def testReadExtendsLeft2(self):
        self.check('''
                                               ..<<<...........
                                             ....>>>...........
                                      ...........>>>....
                   ''')

    def testReadExtendsRight1(self):
        self.check('''
                                             <<<...........
                                             >>>.................
                            .................>>>
                   ''')

    def testReadExtendsRight2(self):
        self.check('''
                                             ....<<<...........
                                             ....>>>.................
                                .................>>>....
                   ''')

    def testReadExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..>>>...
                                     ...>>>..
                  ''')

    def testHitExtendsLeftReadExtendsRight(self):
        self.check('''
                                    ....<<<...........
                                      ..>>>...............
                         ...............>>>..
                  ''')

    def testHitExtendsRightReadExtendsLeft(self):
        self.check('''
                                      ..<<<...............
                                    ....>>>...........
                             ...........>>>....
                  ''')


class HitNegativeReadNegative(TestCase, TestNormalizeHSPMixin):

    """
    This class appears to not be needed. As far as we have seen, the
    read is always positive with ascending start, end offsets.
    """

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

    def testHitExtendsLeft1(self):
        self.check('''
                                             ....<<<
                                               ..<<<
                                               ..>>>
                   ''')

    def testHitExtendsLeft2(self):
        self.check('''
                                             ....<<<...........
                                               ..<<<...........
                                               ..>>>...........
                   ''')

    def testHitExtendsRight1(self):
        self.check('''
                                             <<<.................
                                             <<<...........
                                             <<<...........
                   ''')

    def testHitExtendsRight2(self):
        self.check('''
                                             ....<<<.................
                                             ....<<<...........
                                             ....<<<...........
                   ''')

    def testHitExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..<<<...
                                      ..<<<...
                  ''')

    def testReadExtendsLeft1(self):
        self.check('''
                                               ..<<<
                                             ....<<<
                                             ....<<<
                   ''')

    def testReadExtendsLeft2(self):
        self.check('''
                                               ..<<<...........
                                             ....<<<...........
                                             ....<<<...........
                   ''')

    def testReadExtendsRight1(self):
        self.check('''
                                             <<<...........
                                             <<<.................
                                             <<<.................
                   ''')

    def testReadExtendsRight2(self):
        self.check('''
                                             ....<<<...........
                                             ....<<<.................
                                             ....<<<.................
                   ''')

    def testReadExtendsBoth(self):
        self.check('''
                                  ......<<<...........
                                      ..<<<...
                                      ..<<<...
                  ''')

    def testHitExtendsLeftReadExtendsRight(self):
        self.check('''
                                    ....<<<...........
                                      ..<<<...............
                                      ..<<<...............
                  ''')

    def testHitExtendsRightReadExtendsLeft(self):
        self.check('''
                                      ..<<<...............
                                    ....<<<...........
                                    ....<<<...........
                  ''')


class Old_ReadPositiveHitPositive(TestCase):
    """
    Tests for normalizeHSP when the hit start is less than the hit end.

    NOTE: Please don't add any tests below. Use the Template based tests above.
    """

    frame = Frame(read=1, hit=1)

    def testIdentical(self):
        """The hit start and end are identical to those of the read.

             ssss
             qqqq
        """
        hsp = FakeHSP(hitStart=1, hitEnd=4, readStart=1, readEnd=4,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInHit': 0,
            'readEndInHit': 4,
        }, normalized)

    def testHitExtendsLeft(self):
        """The hit overlaps the read to the left.

              ssssss
                qqqq
        """
        hsp = FakeHSP(hitStart=3, hitEnd=6, readStart=1, readEnd=4,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastn')
        self.assertEqual({
            'hitStart': 2,
            'hitEnd': 6,
            'readStart': 0,
            'readEnd': 4,
            'readStartInHit': 2,
            'readEndInHit': 6,
        }, normalized)

    def testReadExtendsLeft(self):
        """
        The read sticks out to the left of the hit.

               ssss
             qqqqqq
        """
        hsp = FakeHSP(hitStart=1, hitEnd=4, readStart=3, readEnd=6,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInHit': -2,
            'readEndInHit': 4,
        }, normalized)

    def testReadExtendsRight(self):
        """The read sticks out to the right of the hit.

               ssss
               qqqqqq
        """
        hsp = FakeHSP(hitStart=1, hitEnd=4, readStart=1, readEnd=4,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInHit': 0,
            'readEndInHit': 6,
        }, normalized)

    def testReadExtendsRightAndLeft(self):
        """The read extends to the right and left of the hit.

                ssss
               qqqqqq
        """
        hsp = FakeHSP(hitStart=1, hitEnd=4, readStart=2, readEnd=5,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 1,
            'readEnd': 5,
            'readStartInHit': -1,
            'readEndInHit': 5,
        }, normalized)

    def testHitExtendsRightAndLeft(self):
        """The hit extends to the right and left of the read.

               sssssss
                qqqq
        """
        hsp = FakeHSP(hitStart=2, hitEnd=5, readStart=1, readEnd=4,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastn')
        self.assertEqual({
            'hitStart': 1,
            'hitEnd': 5,
            'readStart': 0,
            'readEnd': 4,
            'readStartInHit': 1,
            'readEndInHit': 5,
        }, normalized)

    def test20131115Debugging(self):
        """
        This is an example I manually examined for Barbara on 2013-11-15.
        """

        read = 'TTCTTTTTGCATTTGATAGT-TTGCTACAAG'
        hit = 'TTCTTTTTGCAATAGTCAGTCTTGCTAAAAG'
        hsp = FakeHSP(hitStart=45, hitEnd=75, readStart=120,
                      readEnd=149, frame=self.frame, read=read, hit=hit)
        normalized = normalizeHSP(hsp, 149, 'blastn')
        self.assertEqual({
            'hitStart': 44,
            'hitEnd': 75,
            'readStart': 119,
            'readEnd': 149,
            'readStartInHit': -75,
            'readEndInHit': 75,
        }, normalized)


class Old_ReadPositiveHitNegative(TestCase):
    """
    Tests for normalizeHSP when the hit start is greater than the hit
    end.

    NOTE: Please don't add any tests below. Use the Template based tests above.
    """
    frame = Frame(read=1, hit=-1)

    def testIdentical(self):
        """
        The hit start and end are identical to those of the read.

              ssss
              qqqq
        """
        hsp = FakeHSP(hitStart=4, hitEnd=1, readStart=1, readEnd=4,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInHit': 0,
            'readEndInHit': 4,
        }, normalized)

    def testReadExtendsLeft(self):
        """The read sticks out to the left of the hit.

                ssss
                qqqqqq
        """
        hsp = FakeHSP(hitStart=4, hitEnd=1, readStart=1, readEnd=4,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 0,
            'readEnd': 4,
            'readStartInHit': -2,
            'readEndInHit': 4,
        }, normalized)

    def testReadExtendsLeft2(self):
        """The read sticks out to the left of the hit.

              ssss
                qqqqqq
        """
        hsp = FakeHSP(hitStart=4, hitEnd=1, readStart=3, readEnd=6,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInHit': 0,
            'readEndInHit': 6,
        }, normalized)

    def testReadExtendsRight(self):
        """The read sticks out to the right of the hit.

                  ssss
                qqqqqq
        """
        hsp = FakeHSP(hitStart=4, hitEnd=1, readStart=3, readEnd=6,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInHit': 0,
            'readEndInHit': 6,
        }, normalized)

    def testReadExtendsRight2(self):
        """The read sticks out to the right of the read.

                ssss
            qqqqqq
        """
        hsp = FakeHSP(hitStart=2, hitEnd=1, readStart=5, readEnd=6,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 6, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 2,
            'readStart': 4,
            'readEnd': 6,
            'readStartInHit': 0,
            'readEndInHit': 6,
        }, normalized)

    def testReadExtendsRightAndLeft(self):
        """The read extends to the right and left of the hit.

                  ssss
                qqqqqqq
        """
        hsp = FakeHSP(hitStart=4, hitEnd=1, readStart=3, readEnd=6,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 7, 'blastn')
        self.assertEqual({
            'hitStart': 0,
            'hitEnd': 4,
            'readStart': 2,
            'readEnd': 6,
            'readStartInHit': -1,
            'readEndInHit': 6,
        }, normalized)

    def testHitExtendsRightAndLeft(self):
        """The hit extends to the right and left of the read.

                sssssss
                 qqqq
        """
        hsp = FakeHSP(hitStart=5, hitEnd=2, readStart=1, readEnd=4,
                      frame=self.frame)
        normalized = normalizeHSP(hsp, 4, 'blastn')
        self.assertEqual({
            'hitStart': 1,
            'hitEnd': 5,
            'readStart': 0,
            'readEnd': 4,
            'readStartInHit': 1,
            'readEndInHit': 5,
        }, normalized)

    def test20130721Debugging(self):
        """
        This is an example I manually examined on 2013-07-21.

        I had to invent hit and read strings though on 2013-11-21 to
        introduce 5 gaps in the read due to more rigorous checking in
        normalizeHSP.
        """

        hit = (
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'AAAAA')

        read = (
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGT'
            '-----')

        hsp = FakeHSP(hitStart=9018, hitEnd=8764, readStart=66,
                      readEnd=315, frame=self.frame, hit=hit, read=read)
        normalized = normalizeHSP(hsp, 316, 'blastn')
        self.assertEqual({
            'hitStart': 8763,
            'hitEnd': 9018,
            'readStart': 65,
            'readEnd': 315,
            'readStartInHit': 8762,
            'readEndInHit': 9083,
        }, normalized)

    def test20131113Debugging(self):
        """
        This is an example I manually examined on 2013-11-13.
        """
        hit = (
            'GTCGAGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGTGTCCGGCA'
            'AGGTTGCCAAGGAGCAGATCGACATCGATAACGCCAAGCACACCAAGTGATGCACTGA'
            'CGACGGGTGAGGCCCAGATTCCTACGGCCTGGGCCTCTGTCTGCGTCGGGATGCCATT'
            'AGGCCGGTAGGATCGGTCACATGATCGATCCCAAGCTCCTGCGAACGGATCCGGACGC'
            'CGTTCGTCGCTCCCAGGCCGCCCGCGGCGAGGACTCCTCGGTTGTGGACGACGTTGTC'
            'GCCGCAGATGAGGCTCGTCGTGAGGCTATTGCTGCCCATGAGAACCTGCGTGCAGAAC'
            'AGAAGGGACTCGGCAAGCGAATCGCTAAAGCATCCGGTG')

        read = (
            'GTC-AGAAGATCAAGATTGGTAAGGAGGCCGTGCAGGACACCGAGACCGTGTCCGGCA'
            'AGGTTGCCAAGGAGCAGATCGACATCGATAACGCCAAGCACACCAAGTGATGCACTGA'
            'CGACGGGTGAGGCCCAGATTCCTACGGCCTGGGCCTCTGTCTGCGTCGGGATGCCATT'
            'AGGCCGCTAGGATCGGTCACATGATCGATCCCAAGCTCCTGCGAACGGATCCGGACGC'
            'CGTTCGTCGCTCCCAGGCCGCCCGCGGCGAGGACTCCTCGGTTGTGGACGACGTTGTC'
            'GCCGCAGATGAGGCTCGTCGTGAGGCTATTGCTGCCCATGAGAACCTGCGTGCAGAAC'
            'AGAAGGGACTCGGCAAGCGAATCGCTAAAGCATCCGGTG')

        hsp = FakeHSP(hitStart=2339751, hitEnd=2339365, readStart=1,
                      readEnd=386, frame=self.frame, hit=hit, read=read)
        normalized = normalizeHSP(hsp, 396, 'blastn')
        self.assertEqual({
            'hitStart': 2339364,
            'hitEnd': 2339751,
            'readStart': 0,
            'readEnd': 386,
            'readStartInHit': 2339354,
            'readEndInHit': 2339751,
        }, normalized)

    def test20131115Debugging(self):
        """
        This is an example I manually examined for BM on 2013-11-15.
        """
        read = 'CTCTTGCA-CCTTAGGTACC'
        hit = 'CTCTAGCAGCCTTAGGTACC'
        hsp = FakeHSP(hitStart=1776, hitEnd=1795, readStart=131,
                      readEnd=149, frame=self.frame, read=read, hit=hit)
        normalized = normalizeHSP(hsp, 149, 'blastn')
        self.assertEqual({
            'hitStart': 1775,
            'hitEnd': 1795,
            'readStart': 130,
            'readEnd': 149,
            'readStartInHit': 1775,
            'readEndInHit': 1925,
        }, normalized)

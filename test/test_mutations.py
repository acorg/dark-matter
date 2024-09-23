from unittest import TestCase

from dark.mutations import getAPOBECFrequencies, mutateString


class TestGetAPOBECFrequencies(TestCase):
    """
    Tests the getAPOBECFrequencies function
    """

    def testCorrectResult(self):
        """
        The function should return the right bases.
        s:  AGTCTAAAGTCATGACTGGTCCCCTTG
        q:      ..AA..T....G....G...
        """
        result = []
        subjectTitle = "333031"
        subjectSequence = "AGTCTAAAGTCATGACTGGTCCCCTTG"
        subject = "%s \t \t \t %s" % (subjectTitle, subjectSequence)
        result.append(subject)

        queryTitle = "H9L0NJQ01DXGOQ"
        querySequence = "    ..AA..T....G....G...   "
        query = "%s \t %s" % (queryTitle, querySequence)
        result.append(query)

        bases = getAPOBECFrequencies(result, "C", "G", "cPattern")

        self.assertEqual(
            {
                "ACA": 0,
                "ACC": 0,
                "ACG": 0,
                "ACT": 1,
                "CCA": 0,
                "CCC": 0,
                "CCG": 0,
                "CCT": 0,
                "GCA": 0,
                "GCC": 0,
                "GCG": 0,
                "GCT": 0,
                "TCA": 0,
                "TCC": 1,
                "TCG": 0,
                "TCT": 0,
            },
            bases,
        )


class TestMutateString(TestCase):
    """
    Tests the mutateString function.
    """

    def testZeroLength(self):
        """
        The function should raise when given an empty string
        """
        error = "Empty original string passed."
        with self.assertRaisesRegex(ValueError, error):
            mutateString("", 1)

    def testTooManyMutationsRequested(self):
        """
        The function should raise when asked to perform more mutations than
        there are in the original string.
        """
        error = "Cannot make 2 mutations in a string of length 1"
        with self.assertRaisesRegex(ValueError, error):
            mutateString("x", 2)

    def testDuplicateReplacementLetter(self):
        """
        The function should raise when given a replacement string that contains
        a duplicate letter.
        """
        error = "Replacement string contains duplicates"
        with self.assertRaisesRegex(ValueError, error):
            mutateString("x", 1, "aa")

    def testReplacementLengthOneAppearsInOriginal(self):
        """
        The function should raise when given a replacement string that contains
        just one letter if that letter also appears in the original.
        """
        error = "Impossible replacement"
        with self.assertRaisesRegex(ValueError, error):
            mutateString("x", 1, "x")

    def testZeroReplacements(self):
        """
        The function should return the original string when zero mutations
        are requested.
        """
        self.assertEqual("acgt", mutateString("acgt", 0, "x"))

    def testOneDeterministicReplacement(self):
        """
        The function should return the correct result when one mutation is
        requested and only one mutation is possible.
        """
        self.assertEqual("c", mutateString("a", 1, "ac"))

    def testFiveDeterminsticReplacements(self):
        """
        The function should return the correct result when five mutations are
        requested and only one outcome is possible.
        """
        self.assertEqual("ccccc", mutateString("aaaaa", 5, "ac"))

    def testOneReplacement(self):
        """
        The function should return all possible correct results (and only those
        results) when one mutation is requested and two outcomes are possible.

        NOTE: this test is a bit bad because it's non-deterministic.
        """
        possible = set(["ab", "ba"])
        seen = set()

        for _ in range(100):
            result = mutateString("aa", 1, "ab")
            self.assertTrue(result in possible)
            seen.add(result)
        self.assertEqual(seen, possible)

    def testTwoReplacements(self):
        """
        The function should return all possible correct results (and only those
        results) when two mutations are requested and four outcomes are
        possible.

        NOTE: this test is a bit bad because it's non-deterministic.
        """
        possible = set(["ab", "ba", "ac", "ca"])
        seen = set()

        for _ in range(100):
            result = mutateString("aa", 1, "abc")
            self.assertTrue(result in possible, "%s is not in %r" % (result, possible))
            seen.add(result)
        self.assertEqual(seen, possible)

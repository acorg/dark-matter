from unittest import TestCase

from dark.simplify import simplifyTitle


class SimplifyTitle(TestCase):
    """
    Tests for the dark.simplify.simplifyTitle function.
    """

    def testEmptyTitle(self):
        """
        Simplifying an empty title with a non-empty target should return
        an empty title.
        """
        self.assertEqual('', simplifyTitle('', 'xxx'))

    def testEmtpyTitleWithEmptyTarget(self):
        """
        Simplifying an empty title should return an empty title.
        """
        self.assertEqual('', simplifyTitle('', ''))

    def testPrefix(self):
        """
        When the target is a prefix, the title up to the target (including the
        whole word that has the prefix) should be returned.
        """
        self.assertEqual(
            'Funny sea lion polyoma',
            simplifyTitle('Funny sea lion polyomavirus 1 CSL6994', 'polyoma'))

    def testSuffix(self):
        """
        When the target is a suffix, the title up to the target (including the
        whole word that has the suffix) should be returned.
        """
        self.assertEqual(
            'Funny sea lion polyomavirus',
            simplifyTitle('Funny sea lion polyomavirus 1 CSL6994', 'virus'))

    def testContained(self):
        """
        When the target is contained, the title up to the target (including the
        prefix of the word that has the target) should be returned.
        """
        self.assertEqual(
            'Funny sea lion polyoma',
            simplifyTitle('Funny sea lion polyomavirus 1 CSL6994', 'yoma'))

    def testExact(self):
        """
        When the target is the same as a word in the title, the title up to
        and including the target should be returned.
        """
        self.assertEqual(
            'Funny sea lion',
            simplifyTitle('Funny sea lion polyomavirus 1 CSL6994', 'lion'))

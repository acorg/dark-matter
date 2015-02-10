from unittest import TestCase

from dark.taxonomy import LineageFetcher


class FakeCursor(object):
    def __init__(self, results):
        self._results = results
        self._index = -1

    def execute(self, p):
        pass

    def fetchone(self):
        self._index += 1
        return self._results[self._index]

    def close(self):
        pass


class FakeDbConnection(object):
    def __init__(self, results):
        self._results = results
        self.open = True

    def cursor(self):
        return FakeCursor(self._results)

    def close(self):
        self.open = False


class TestLineageFetcher(TestCase):
    """
    Test LineageFetcher class.
    """
    def testGetTaxonomy(self):
        """
        Test if the LineageFetcher class works properly.
        """
        title = 'gi|5|gb|EU375804.1| Merkel cell polyomavirus'

        db = FakeDbConnection([[15], [2], ['Merkel cell polyomavirus'],
                               [3], ['Polyomavirus'], [2],
                               ['dsDNA viruses'], [1], ['Vira']])
        cursor = db.cursor()

        lineageFetcher = LineageFetcher(db=db, cursor=cursor)

        lineage = lineageFetcher.lineage(title)
        self.assertEqual(['Merkel cell polyomavirus', 'Polyomavirus',
                          'dsDNA viruses', 'Vira'], lineage)

from unittest import TestCase
import six
import sqlite3

from dark.taxonomy import (
    LineageFetcher, Taxonomy, isRetrovirus, isRNAVirus)


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

        db = FakeDbConnection([
            [15], ['Merkel cell polyomavirus'],
            [4], ['Polyomavirus'],
            [3], ['dsDNA viruses'],
            [2], ['Vira'],
            [1],
        ])
        cursor = db.cursor()

        lineageFetcher = LineageFetcher(db=db, cursor=cursor)

        lineage = lineageFetcher.lineage(title)
        self.assertEqual(
            [
                (15, 'Merkel cell polyomavirus'),
                (4, 'Polyomavirus'),
                (3, 'dsDNA viruses'),
                (2, 'Vira'),
            ],
            lineage)


class TestTaxonomy(TestCase):
    """
    Test the Taxonomy class.
    """

    def _makeFetcher(self, data):
        """
        Make a taxonomy database and put it into an Taxonomy
        instance.

        @param data: A C{dict} with 'taxids', 'names', and 'nodes' keys,
            each of which is a C{list} of C{tuple}s containing values to
            insert into the respective tables.
        @return: A C{Taxonomy} instance.
        """
        db = sqlite3.connect(':memory:')
        cursor = db.cursor()

        # The tables created here must be identical to those expected by
        # the Taxonomy class. See also
        # https://github.com/acorg/ncbi-taxonomy-database
        cursor.executescript('''
            CREATE TABLE nodes (
                taxid INTEGER NOT NULL,
                parent_taxid INTEGER NOT NULL,
                rank VARCHAR NOT NULL
            );

            CREATE TABLE accession_taxid (
                accession VARCHAR UNIQUE PRIMARY KEY,
                taxid INTEGER NOT NULL
            );

            CREATE TABLE names (
                taxid INTEGER NOT NULL,
                name VARCHAR NOT NULL
            );

            CREATE INDEX nodes_idx ON nodes(taxid);
            CREATE INDEX accession_idx ON accession_taxid(accession);
            CREATE INDEX name_idx ON names(taxid);
        ''')

        execute = cursor.execute

        for accession, taxid in data.get('taxids', []):
            execute('INSERT INTO accession_taxid VALUES (?, ?)',
                    (accession, taxid))

        for taxid, name in data.get('names', []):
            execute('INSERT INTO names VALUES (?, ?)',
                    (taxid, name))

        for taxid, parentTaxid, rank in data.get('nodes', []):
            execute('INSERT INTO nodes VALUES (?, ?, ?)',
                    (taxid, parentTaxid, rank))

        return Taxonomy(db)

    def testUnknownTaxid(self):
        """
        If a taxonomy id is not present in the accession_taxid table, the
        lineage fetcher must return C{None}.
        """
        fetcher = self._makeFetcher({})
        self.assertIsNone(fetcher.lineage('DQ011818.1'))
        fetcher.close()

    def testReusingClosedFetcher(self):
        """
        Trying to re-use a closed fetcher must result in an AttributeError
        being raised.
        """
        fetcher = self._makeFetcher({})
        fetcher.close()
        error = "^'NoneType' object has no attribute 'cursor'$"
        six.assertRaisesRegex(self, AttributeError, error, fetcher.lineage,
                              'DQ011818.1')

    def testTaxidNotInNamesTable(self):
        """
        If a taxonomy id is not present in the names table, a ValueError must
        be raised.
        """
        fetcher = self._makeFetcher({
            'taxids': (
                ('DQ011818.1', 500),
            ),
        })
        error = '^Could not find taxonomy id 500 in names table$'
        six.assertRaisesRegex(self, ValueError, error, fetcher.lineage,
                              'DQ011818.1')
        fetcher.close()

    def testTaxidNotInNodesTable(self):
        """
        If a taxonomy id is not present in the nodes table, a ValueError must
        be raised.
        """
        fetcher = self._makeFetcher({
            'taxids': (
                ('DQ011818.1', 500),
            ),
            'names': (
                (500, 'Porcine endogenous retrovirus A'),
            ),
        })
        error = '^Could not find taxonomy id 500 in nodes table$'
        six.assertRaisesRegex(self, ValueError, error, fetcher.lineage,
                              'DQ011818.1')
        fetcher.close()

    def testFullLookup(self):
        """
        The full taxonomy of an accession number must be retrievable.
        """
        fetcher = self._makeFetcher({
            'taxids': (
                ('DQ011818.1', 500),
            ),
            'names': (
                (500, 'Porcine endogenous retrovirus A'),
                (400, 'Retroviruses'),
                (300, 'Viruses'),
            ),
            'nodes': (
                (500, 400, 'species'),
                (400, 300, 'genus'),
                (300, 1, 'realm')
            ),
        })

        self.assertEqual(
            (
                (500, 'Porcine endogenous retrovirus A', 'species'),
                (400, 'Retroviruses', 'genus'),
                (300, 'Viruses', 'realm'),
            ),
            fetcher.lineage('DQ011818.1'))
        fetcher.close()


class TestIsRetrovirus(TestCase):
    """
    Test the isRetrovirus function.
    """
    def testYes(self):
        """
        A retrovirus must be detected.
        """
        self.assertTrue(
            isRetrovirus((
                (11676, 'Human immunodeficiency virus 1', 'species'),
                (11646, 'Lentivirus', 'genus'),
                (327045, 'Orthoretrovirinae', 'subfamily'),
                (11632, 'Retroviridae', 'family'),
                (2169561, 'Ortervirales', 'order'),
                (10239, 'Viruses', 'superkingdom'),
            )))

    def testNo(self):
        """
        A non-retrovirus must be detected.
        """
        self.assertFalse(
            isRetrovirus((
                (10310, 'Human alphaherpesvirus 2', 'species'),
                (10294, 'Simplexvirus', 'genus'),
                (10293, 'Alphaherpesvirinae', 'subfamily'),
                (10292, 'Herpesviridae', 'family'),
                (548681, 'Herpesvirales', 'order'),
                (10239, 'Viruses', 'superkingdom'),
            )))


class TestIsRNAVirus(TestCase):
    """
    Test the isRNAVirus function.
    """
    def testByRealm(self):
        """
        Anything in the Riboviria realm must be classified as an RNA virus.
        """
        self.assertTrue(
            isRNAVirus((
                (11234, 'Measles morbillivirus', 'species'),
                (11229, 'Morbillivirus', 'genus'),
                (2560076, 'Orthoparamyxovirinae', 'subfamily'),
                (11158, 'Paramyxoviridae', 'family'),
                (11157, 'Mononegavirales', 'order'),
                (2497574, 'Monjiviricetes', 'class'),
                (2497570, 'Haploviricotina', 'subphylum'),
                (2497569, 'Negarnaviricota', 'phylum'),
                (2559587, 'Riboviria', 'realm'),
                (10239, 'Viruses', 'superkingdom'),
            )))

    def testHBV(self):
        """
        An HBV virus must not be classified as an RNA virus.
        """
        self.assertFalse(
            isRNAVirus((
                (1508712, 'Tent-making bat hepatitis B virus', 'species'),
                (10405, 'Orthohepadnavirus', 'genus'),
                (10404, 'Hepadnaviridae', 'family'),
                (10239, 'Viruses', 'superkingdom'),
            )))

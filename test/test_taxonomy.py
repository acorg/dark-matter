from unittest import TestCase
import sqlite3

from dark.taxonomy import (
    LineageElement as LE,
    LineageFetcher,
    Taxonomy,
    isRetrovirus,
    isRNAVirus,
    isAllowedTaxonomicRank,
)


class FakeCursor:
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


class FakeDbConnection:
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
        title = "gi|5|gb|EU375804.1| Merkel cell polyomavirus"

        db = FakeDbConnection(
            [
                [15],
                ["Merkel cell polyomavirus"],
                [4],
                ["Polyomavirus"],
                [3],
                ["dsDNA viruses"],
                [2],
                ["Vira"],
                [1],
            ]
        )
        cursor = db.cursor()

        lineageFetcher = LineageFetcher(db=db, cursor=cursor)

        lineage = lineageFetcher.lineage(title)
        self.assertEqual(
            [
                (15, "Merkel cell polyomavirus"),
                (4, "Polyomavirus"),
                (3, "dsDNA viruses"),
                (2, "Vira"),
            ],
            lineage,
        )


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
        db = sqlite3.connect(":memory:")
        cursor = db.cursor()

        # The tables created here must be identical to those expected by
        # the Taxonomy class. See also
        # https://github.com/acorg/ncbi-taxonomy-database
        cursor.executescript(
            """
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
        """
        )

        execute = cursor.execute

        for accession, taxid in data.get("taxids", []):
            execute("INSERT INTO accession_taxid VALUES (?, ?)", (accession, taxid))

        for taxid, name in data.get("names", []):
            execute("INSERT INTO names VALUES (?, ?)", (taxid, name))

        for taxid, parentTaxid, rank in data.get("nodes", []):
            execute("INSERT INTO nodes VALUES (?, ?, ?)", (taxid, parentTaxid, rank))

        return Taxonomy(db)

    def testUnknownTaxid(self):
        """
        If a taxonomy id is not present in the accession_taxid table, the
        lineage fetcher must raise a C{ValueError}.
        """
        fetcher = self._makeFetcher({})
        error = r"^Could not find 'DQ011818\.1' in accession_taxid or names " r"tables$"
        self.assertRaisesRegex(ValueError, error, fetcher.lineage, "DQ011818.1")
        fetcher.close()

    def testReusingClosedFetcher(self):
        """
        Trying to re-use a closed fetcher must result in an AttributeError
        being raised.
        """
        fetcher = self._makeFetcher({})
        fetcher.close()
        error = "^'NoneType' object has no attribute 'cursor'$"
        self.assertRaisesRegex(AttributeError, error, fetcher.lineage, "DQ011818.1")

    def testTaxidNotInNamesTable(self):
        """
        If a taxonomy id is not present in the names table, a ValueError must
        be raised.
        """
        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
            }
        )
        error = "^Could not find taxonomy id 500 in names table$"
        self.assertRaisesRegex(ValueError, error, fetcher.lineage, "DQ011818.1")
        fetcher.close()

    def testTaxidNotInNodesTable(self):
        """
        If a taxonomy id is not present in the nodes table, a ValueError must
        be raised.
        """
        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
                "names": ((500, "Porcine endogenous retrovirus A"),),
            }
        )
        error = "^Could not find taxonomy id 500 in nodes table$"
        self.assertRaisesRegex(ValueError, error, fetcher.lineage, "DQ011818.1")
        fetcher.close()

    def testLookupByAccession(self):
        """
        The full taxonomy of an accession number must be retrievable.
        """
        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
                "names": (
                    (500, "Porcine endogenous retrovirus A"),
                    (400, "Retroviruses"),
                    (300, "Viruses"),
                ),
                "nodes": (
                    (500, 400, "species"),
                    (400, 300, "genus"),
                    (300, 1, "realm"),
                ),
            }
        )

        self.assertEqual(
            (
                (500, "Porcine endogenous retrovirus A", "species"),
                (400, "Retroviruses", "genus"),
                (300, "Viruses", "realm"),
            ),
            fetcher.lineage("DQ011818.1"),
        )
        fetcher.close()

    def testLookupByName(self):
        """
        The taxonomy must be retrievable by name.
        """
        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
                "names": (
                    (500, "Porcine endogenous retrovirus A"),
                    (400, "Retroviruses"),
                    (300, "Viruses"),
                ),
                "nodes": (
                    (500, 400, "species"),
                    (400, 300, "genus"),
                    (300, 1, "realm"),
                ),
            }
        )

        self.assertEqual(
            (
                (500, "Porcine endogenous retrovirus A", "species"),
                (400, "Retroviruses", "genus"),
                (300, "Viruses", "realm"),
            ),
            fetcher.lineage("Porcine endogenous retrovirus A"),
        )
        fetcher.close()

    def testLookupByTaxid(self):
        """
        The taxonomy must be retrievable by taxonomy id.
        """
        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
                "names": (
                    (500, "Porcine endogenous retrovirus A"),
                    (400, "Retroviruses"),
                    (300, "Viruses"),
                ),
                "nodes": (
                    (500, 400, "species"),
                    (400, 300, "genus"),
                    (300, 1, "realm"),
                ),
            }
        )

        self.assertEqual(
            (
                (500, "Porcine endogenous retrovirus A", "species"),
                (400, "Retroviruses", "genus"),
                (300, "Viruses", "realm"),
            ),
            fetcher.lineage(500),
        )
        fetcher.close()

    def testLookupWithSkip(self):
        """
        It must be possible to skip a level when retrieving a taxonomy.
        """

        def skipFunc(lineageElement):
            return lineageElement.rank == "genus"

        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
                "names": (
                    (500, "Porcine endogenous retrovirus A"),
                    (400, "Retroviruses"),
                    (300, "Viruses"),
                ),
                "nodes": (
                    (500, 400, "species"),
                    (400, 300, "genus"),
                    (300, 1, "realm"),
                ),
            }
        )

        self.assertEqual(
            (
                (500, "Porcine endogenous retrovirus A", "species"),
                (300, "Viruses", "realm"),
            ),
            fetcher.lineage(500, skipFunc=skipFunc),
        )
        fetcher.close()

    def testLookupWithStop(self):
        """
        It must be possible to stop at a level when retrieving a taxonomy.
        """

        def stopFunc(lineageElement):
            return lineageElement.rank == "genus"

        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
                "names": (
                    (500, "Porcine endogenous retrovirus A"),
                    (400, "Retroviruses"),
                    (300, "Viruses"),
                ),
                "nodes": (
                    (500, 400, "species"),
                    (400, 300, "genus"),
                    (300, 1, "realm"),
                ),
            }
        )

        self.assertEqual(
            (
                (500, "Porcine endogenous retrovirus A", "species"),
                (400, "Retroviruses", "genus"),
            ),
            fetcher.lineage(500, stopFunc=stopFunc),
        )
        fetcher.close()

    def testLookupWithSkipAndStop(self):
        """
        It must be possible to stop at a level when retrieving a taxonomy
        and skip the last (stopping) level.
        """

        def skipFunc(lineageElement):
            return lineageElement.rank == "genus"

        def stopFunc(lineageElement):
            return lineageElement.rank == "genus"

        fetcher = self._makeFetcher(
            {
                "taxids": (("DQ011818.1", 500),),
                "names": (
                    (500, "Porcine endogenous retrovirus A"),
                    (400, "Retroviruses"),
                    (300, "Viruses"),
                ),
                "nodes": (
                    (500, 400, "species"),
                    (400, 300, "genus"),
                    (300, 1, "realm"),
                ),
            }
        )

        self.assertEqual(
            ((500, "Porcine endogenous retrovirus A", "species"),),
            fetcher.lineage(500, skipFunc=skipFunc, stopFunc=stopFunc),
        )
        fetcher.close()

    def testSubsetLineageByRanks(self):
        """
        Test the Taxonomy.subsetLineageByRanks static method.
        """
        lineage = (
            LE(11234, "Measles morbillivirus", "species"),
            LE(11229, "Morbillivirus", "genus"),
            LE(2560076, "Orthoparamyxovirinae", "subfamily"),
            LE(11158, "Paramyxoviridae", "family"),
            LE(11157, "Mononegavirales", "order"),
            LE(2497574, "Monjiviricetes", "class"),
            LE(2497570, "Haploviricotina", "subphylum"),
            LE(2497569, "Negarnaviricota", "phylum"),
            LE(2559587, "Riboviria", "realm"),
            LE(10239, "Viruses", "superkingdom"),
        )

        def filterFunc(rank):
            return rank in {"phylum", "genus"}

        self.assertEqual(
            (
                LE(11229, "Morbillivirus", "genus"),
                LE(2497569, "Negarnaviricota", "phylum"),
            ),
            tuple(Taxonomy.subsetLineageByRanks(lineage, filterFunc)),
        )


class TestIsRetrovirus(TestCase):
    """
    Test the isRetrovirus function.
    """

    def testYes(self):
        """
        A retrovirus must be detected.
        """
        self.assertTrue(
            isRetrovirus(
                (
                    LE(11676, "Human immunodeficiency virus 1", "species"),
                    LE(11646, "Lentivirus", "genus"),
                    LE(327045, "Orthoretrovirinae", "subfamily"),
                    LE(11632, "Retroviridae", "family"),
                    LE(2169561, "Ortervirales", "order"),
                    LE(10239, "Viruses", "superkingdom"),
                )
            )
        )

    def testNo(self):
        """
        A non-retrovirus must be detected.
        """
        self.assertFalse(
            isRetrovirus(
                (
                    LE(10310, "Human alphaherpesvirus 2", "species"),
                    LE(10294, "Simplexvirus", "genus"),
                    LE(10293, "Alphaherpesvirinae", "subfamily"),
                    LE(10292, "Herpesviridae", "family"),
                    LE(548681, "Herpesvirales", "order"),
                    LE(10239, "Viruses", "superkingdom"),
                )
            )
        )


class TestIsRNAVirus(TestCase):
    """
    Test the isRNAVirus function.
    """

    def testByRealm(self):
        """
        Anything in the Riboviria realm must be classified as an RNA virus.
        """
        self.assertTrue(
            isRNAVirus(
                (
                    LE(11234, "Measles morbillivirus", "species"),
                    LE(11229, "Morbillivirus", "genus"),
                    LE(2560076, "Orthoparamyxovirinae", "subfamily"),
                    LE(11158, "Paramyxoviridae", "family"),
                    LE(11157, "Mononegavirales", "order"),
                    LE(2497574, "Monjiviricetes", "class"),
                    LE(2497570, "Haploviricotina", "subphylum"),
                    LE(2497569, "Negarnaviricota", "phylum"),
                    LE(2559587, "Riboviria", "realm"),
                    LE(10239, "Viruses", "superkingdom"),
                )
            )
        )

    def testHBV(self):
        """
        An HBV virus must not be classified as an RNA virus.
        """
        self.assertFalse(
            isRNAVirus(
                (
                    LE(1508712, "Tent-making bat hepatitis B virus", "species"),
                    LE(10405, "Orthohepadnavirus", "genus"),
                    LE(10404, "Hepadnaviridae", "family"),
                    LE(10239, "Viruses", "superkingdom"),
                )
            )
        )

    def testHIV(self):
        """
        An HIV virus must be classified as an RNA virus.
        """
        self.assertTrue(
            isRNAVirus(
                (
                    LE(11676, "Human immunodeficiency virus 1", "species"),
                    LE(11646, "Lentivirus", "genus"),
                    LE(327045, "Orthoretrovirinae", "subfamily"),
                    LE(11632, "Retroviridae", "family"),
                    LE(2169561, "Ortervirales", "order"),
                    LE(10239, "Viruses", "superkingdom"),
                )
            )
        )


class TestAllowedTaxonomicRank(TestCase):
    """
    Test the isAllowedTaxonomicRank function.
    """

    def testNoRanksAllowed(self):
        """
        The function must return False when no taxonomic ranks are allowed.
        """
        allowed = set()
        lineage = (
            LE(245, "x", "realm"),
            LE(246, "bacteria", "kingdom"),
            LE(247, "y", "species"),
        )
        self.assertFalse(isAllowedTaxonomicRank(allowed, lineage))

    def testEmptyLineage(self):
        """
        The function must return False when the passed lineage is empty.
        """
        allowed = set((("bacteria", "kingdom"),))
        lineage = ()
        self.assertFalse(isAllowedTaxonomicRank(allowed, lineage))

    def testExactCase(self):
        """
        The function must return True when case matches exactly.
        """
        allowed = set((("bacteria", "kingdom"),))
        lineage = (
            LE(245, "x", "realm"),
            LE(246, "bacteria", "kingdom"),
            LE(247, "y", "species"),
        )
        self.assertTrue(isAllowedTaxonomicRank(allowed, lineage))

    def testNonMatchingCase(self):
        """
        The function must return True when case doesn't match.
        """
        allowed = set((("BACTERIA", "KINGDOM"),))
        lineage = (
            LE(245, "x", "realm"),
            LE(246, "bacteria", "kingdom"),
            LE(247, "y", "species"),
        )
        self.assertTrue(isAllowedTaxonomicRank(allowed, lineage))

    def testMultipleMatches(self):
        """
        The function must return True when more than one part of the lineage matches.
        """
        allowed = set((("bacteria", "kingdom"), ("y", "species")))
        lineage = (
            LE(245, "x", "realm"),
            LE(246, "bacteria", "kingdom"),
            LE(247, "y", "species"),
        )
        self.assertTrue(isAllowedTaxonomicRank(allowed, lineage))

    def testRankMismatch(self):
        """
        The function must return False when the name matches but the rank does not.
        """
        allowed = set((("bacteria", "species"),))
        lineage = (
            LE(245, "x", "realm"),
            LE(246, "bacteria", "kingdom"),
            LE(247, "y", "species"),
        )
        self.assertFalse(isAllowedTaxonomicRank(allowed, lineage))

    def testNameMismatch(self):
        """
        The function must return False when the rank matches but the name does not.
        """
        allowed = set((("virus", "kingdom"),))
        lineage = (
            LE(245, "x", "realm"),
            LE(246, "bacteria", "kingdom"),
            LE(247, "y", "species"),
        )
        self.assertFalse(isAllowedTaxonomicRank(allowed, lineage))

    def testNameAndRankMismatch(self):
        """
        The function must return False when the name and rank both do not match.
        """
        allowed = set((("virus", "genus"),))
        lineage = (
            LE(245, "x", "realm"),
            LE(246, "bacteria", "kingdom"),
            LE(247, "y", "species"),
        )
        self.assertFalse(isAllowedTaxonomicRank(allowed, lineage))


class TestLineageElement(TestCase):
    """
    Test the LineageElement named tuple.
    """

    def testTaxid(self):
        """
        The taxid attribute must be set as expected.
        """
        element = LE(245, "no name", "species")
        self.assertEqual(245, element.taxid)
        self.assertEqual(245, element[0])

    def testName(self):
        """
        The name attribute must be set as expected.
        """
        element = LE(245, "no name", "species")
        self.assertEqual("no name", element.name)
        self.assertEqual("no name", element[1])

    def testRank(self):
        """
        The rank attribute must be set as expected.
        """
        element = LE(245, "no name", "species")
        self.assertEqual("species", element.rank)
        self.assertEqual("species", element[2])

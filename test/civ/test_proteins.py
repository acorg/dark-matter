from unittest import TestCase

# from six import assertRaisesRegex

from dark.civ.proteins import SqliteIndex, SqliteIndexWriter


class TestSqliteIndex(TestCase):
    """
    Test the SqliteIndex class.
    """

    def testCountProteinsEmpty(self):
        """
        An empty database must have zero proteins.
        """
        writer = SqliteIndexWriter(":memory:")
        db = SqliteIndex(writer._connection)
        self.assertEqual(0, db.proteinCount())
        writer.close()
        db.close()

    def testCountGenomesEmpty(self):
        """
        An empty database must have zero genomes.
        """
        writer = SqliteIndexWriter(":memory:")
        db = SqliteIndex(writer._connection)
        self.assertEqual(0, db.genomeCount())
        writer.close()
        db.close()

    def testCountProteinsOne(self):
        """
        A database with one protein must have a protein count of one.
        """
        writer = SqliteIndexWriter(":memory:")
        writer.addProtein("NC222", "NC23", "AAA", "offets", True, False, 3)
        db = SqliteIndex(writer._connection)
        self.assertEqual(1, db.proteinCount())
        writer.close()
        db.close()

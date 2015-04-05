from dark import mysql


class LineageFetcher(object):
    """
    Provide access to the NCBI taxonomy database so we can retrieve the lineage
    of title sequences hit by BLAST.
    """
    def __init__(self, db=None, cursor=None):
        self._db = db or mysql.getDatabaseConnection()
        self._cursor = cursor or self._db.cursor()
        self._cache = {}

    def lineage(self, title):
        """
        For a give title, gets the lineage information from taxonomy
        database.

        @param title: the C{str} title of a sequence hit by BLAST.

        @return: A C{list} of the taxonomic categories of the title. If
            no taxonomy is found, the list will be empty.
        """
        if title in self._cache:
            return self._cache[title]

        lineage = []
        gi = int(title.split('|')[1])
        query = 'SELECT taxID from gi_taxid where gi = %d' % gi
        self._cursor.execute(query)

        try:
            taxID = self._cursor.fetchone()[0]
            while taxID != 1:
                query = ('SELECT parent_taxID from nodes where taxID = %s' %
                         taxID)
                self._cursor.execute(query)
                parentTaxID = self._cursor.fetchone()[0]
                query = 'SELECT name from names where taxId = %s' % taxID
                self._cursor.execute(query)
                scientificName = self._cursor.fetchone()[0]
                lineage.append(scientificName)
                taxID = parentTaxID
        except TypeError:
            lineage = []

        self._cache[title] = lineage
        return lineage

    def close(self):
        """
        Close the database connection and render self invalid. Any subsequent
        re-use of self will raise an error.
        """
        self._cursor.close()
        self._db.close()
        self._cursor = self._db = self._cache = None

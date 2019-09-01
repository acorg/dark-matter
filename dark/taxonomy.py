import sqlite3
from six import string_types
from operator import attrgetter
from cachetools import LRUCache, cachedmethod

from dark.database import getDatabaseConnection


class LineageFetcher(object):
    """
    Provide access to the NCBI taxonomy database so we can retrieve the lineage
    of title sequences hit by BLAST.
    """
    def __init__(self, db=None, cursor=None):
        self._db = db or getDatabaseConnection()
        self._cursor = cursor or self._db.cursor()
        self._cache = {}

    def lineage(self, title):
        """
        Get lineage information from the taxonomy database for a given title.

        @param title: A C{str} sequence title (e.g., from a BLAST hit). Of the
            form 'gi|63148399|gb|DQ011818.1| Description...'. It is the gi
            number (63148399 in this example) that is looked up in the taxonomy
            database.
        @return: A C{list} of the taxonomic categories of the title. Each list
            element is an (C{int}, C{str}) 2-tuple, giving a taxonomy id and
            a scientific name. The first element in the list will correspond to
            C{title}, and each successive element is the parent of the
            preceeding one. If no taxonomy is found, the returned list will be
            empty.
        """
        if title in self._cache:
            return self._cache[title]

        lineage = []
        gi = int(title.split('|')[1])
        query = 'SELECT taxID FROM gi_taxid WHERE gi = %d' % gi

        try:
            while True:
                self._cursor.execute(query)
                taxID = self._cursor.fetchone()[0]
                if taxID == 1:
                    break
                query = 'SELECT name FROM names WHERE taxId = %s' % taxID
                self._cursor.execute(query)
                scientificName = self._cursor.fetchone()[0]
                lineage.append((taxID, scientificName))
                # Move up to the parent.
                query = ('SELECT parent_taxID FROM nodes WHERE taxID = %s' %
                         taxID)
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


class AccessionLineageFetcher(object):
    """
    Provide access to the NCBI taxonomy database so we can retrieve the
    taxonomy lineage corresponding to an accession number.

    @param dbFilenameOrConnection: Either a C{str} database filename or an
        open sqlite3 database. The database must contain taxonomy information
        with tables and columns named as in the building scripts used in
        https://github.com/acorg/ncbi-taxonomy-database
    """
    CACHE_SIZE = 10E6

    def __init__(self, dbFilenameOrConnection):
        if isinstance(dbFilenameOrConnection, string_types):
            self._db = sqlite3.connect(dbFilenameOrConnection)
            self._closeConnection = True
        else:
            self._db = dbFilenameOrConnection
            self._closeConnection = False
        self._cache = LRUCache(maxsize=self.CACHE_SIZE)

    def lineageFromTaxid(self, taxid):
        cursor = self._db.cursor()
        execute = cursor.execute
        fetchone = cursor.fetchone

        lineage = []

        while taxid != 1:
            execute('SELECT name FROM names WHERE taxid = ?', (taxid,))
            result = fetchone()
            if result is None:
                raise ValueError(
                    'Could not find taxonomy id %r in names table' %
                    (taxid,))
            name = result[0]
            execute('SELECT parent_taxid, rank FROM nodes WHERE taxid = ?',
                    (taxid,))
            result = fetchone()
            if result is None:
                raise ValueError(
                    'Could not find taxonomy id %r in nodes table' %
                    (taxid,))
            parentTaxid, rank = result
            lineage.append((taxid, name, rank))
            taxid = parentTaxid

        return tuple(lineage) or None

    @cachedmethod(attrgetter('_cache'))
    def lineage(self, accession):
        """
        Get lineage information from the taxonomy database for an
        accession number or taxonomy id.

        @param accession: A C{str} accession number or C{int} taxonomy id.
            If a C{str}, C{accession} must include the version (e.g., the
            ".1" part of "DQ011818.1") if the taxonomy database in use uses
            it (as is the case with the database built by
            https://github.com/acorg/ncbi-taxonomy-database
        @raise ValueError: If a taxonomy id cannot be found in the names or
            nodes table.
        @raise AttributeError: If C{self.close} has been called.
        @return: A C{tuple} of the taxonomic categories of the title. Each
            tuple element is a 3-tuple of (C{int}, C{str}, C{str}) giving a
            taxonomy id a (scientific) name, and the rank (species, genus,
            etc). The first element in the list corresponds to the passed
            C{accession}, number and each successive element is the parent of
            the preceeding one. The top level of the taxonomy hierarchy
            (with taxid = 1) is not included in the returned list. If no
            taxonomy information is found for C{accession}, return C{None}.
        """
        if isinstance(accession, int):
            return self.lineageFromTaxid(accession)
        else:
            cursor = self._db.cursor()
            fetchone = cursor.fetchone

            cursor.execute(
                'SELECT taxid FROM accession_taxid WHERE accession = ?',
                (accession,))
            row = fetchone()
            if row:
                return self.lineageFromTaxid(int(row[0]))

    def close(self):
        """
        Close the database connection (if we opened it).
        """
        if self._closeConnection:
            self._db.close()
        # Set self._db to None so self.lineage will raise an AttributeError
        # exception if it is called again.
        self._db = None

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, traceback):
        self.close()


def isRetrovirus(lineage):
    """
    Determine whether a lineage corresponds to a retrovirus.

    @param lineage: A C{tuple} of taxonomy id, scientific name, and rank
        as returned by C{AccessionLineageFetcher.lineage}.
    @return: C{True} if the lineage corresponds to a retrovirus, C{False}
        otherwise.
    """
    for taxid, name, rank in lineage:
        if rank == 'family' and name == 'Retroviridae':
            return True

    return False


def isRNAVirus(lineage):
    """
    Determine whether a lineage corresponds to an RNA virus.

    @param lineage: A C{tuple} of taxonomy id, scientific name, and rank
        as returned by C{AccessionLineageFetcher.lineage}.
    @return: C{True} if the lineage corresponds to an RNA virus, C{False}
        otherwise.
    """
    for taxid, name, rank in lineage:
        if rank == 'realm' and name == 'Riboviria':
            return True

    return False


def _preprocessLineage(lineage):
    """
    Pre-process a lineage to make it easier to deal with by others.

    @param lineage: A C{tuple} of taxonomy id, scientific name, and rank
        as returned by C{AccessionLineageFetcher.lineage}.
    @return: A 3-tuple of taxids, names, and ranks (all as C{str}ings).
    """
    taxids, names, ranks = [], [], []

    for (taxid, name, rank) in lineage:
        taxids.append(str(taxid))
        ranks.append('-' if rank == 'no rank' else rank)
        names.append(name)

    return taxids, names, ranks


def formatLineage(lineage, namesOnly=False, separator=None, prefix=''):
    """
    Format a lineage for printing.

    @param lineage: A C{tuple} of taxonomy id, scientific name, and rank
        as returned by C{AccessionLineageFetcher.lineage}.
    @param namesOnly: If C{True} only print taxonomic names.
    @param separator: A C{str} separator to put between fields. If C{None},
        return a space-padded aligned columns.
    @param prefix: A C{str} to put at the start of each line.
    @return: A formatted C{str} for printing.
    """
    taxids, names, ranks = _preprocessLineage(lineage)

    if namesOnly:
        # The separator is guaranteed to be set by our caller.
        return prefix + separator.join(names)

    if separator is not None:
        return '\n'.join(
            '%s%s%s%s%s%s' % (prefix, rank, separator, name, separator, taxid)
            for (taxid, name, rank) in zip(taxids, names, ranks))
    else:
        taxidWidth = max(len(x) for x in taxids)
        nameWidth = max(len(x) for x in names)
        rankWidth = max(len(x) for x in ranks)

        return '\n'.join(
            '%s%-*s %-*s %*s' % (
                prefix, rankWidth, rank, nameWidth, name, taxidWidth, taxid)
            for (taxid, name, rank) in zip(taxids, names, ranks))


def lineageTaxonomyLinks(lineage):
    """
    Get HTML links for a lineage.

    @param lineage: A C{tuple} of taxonomy id, scientific name, and rank
        as returned by C{AccessionLineageFetcher.lineage}.
    @return: A C{list} of HTML C{str} links.
    """
    taxids, names, _ = _preprocessLineage(lineage)

    assert names[-1] == 'Viruses'
    names[0] = 'taxon'

    return [
        '<a href="%s%s">%s</a>' % (
            'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=',
            taxid, name) for taxid, name in list(zip(taxids, names))[:-1]
    ]

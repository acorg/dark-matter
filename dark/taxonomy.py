from __future__ import annotations

from operator import attrgetter
from cachetools import LRUCache, cachedmethod
from json import dumps
import re
from collections import namedtuple
from os import environ
from typing import Optional, Iterable
import argparse

from dark.database import getDatabaseConnection
from dark.sqlite3 import sqliteConnect

TAXONOMY_DATABASE_ENV_VAR = "DARK_MATTER_TAXONOMY_DATABASE"
TAXONOMY_DATABASE_COMMAND_LINE_OPTION = "taxonomyDatabase"

LineageElement = namedtuple("LineageElement", ("taxid", "name", "rank"))

# These are name/rank tuples that indicate that a lineage is from an RNA
# virus.
RNA_VIRUS_LINEAGE_ELEMENTS = set(
    (
        ("Retroviridae", "family"),
        ("Pseudoviridae", "family"),
        ("Riboviria", "realm"),
    )
)

FUNGUS_ONLY_VIRUS_REGEX = re.compile(r"\b(?:mycovirus)\b", re.I)

FUNGUS_ONLY_GENERA: set[str] = set()

FUNGUS_ONLY_FAMILIES = {
    "Chrysoviridae",  # Infects penicillum. E.g., NC_040738.1
    "Totiviridae",  # Actually infects fungi and protozoa.
}

PLANT_ONLY_VIRUS_REGEX = re.compile(
    r"\b(?:grapevine|mosaic|mottle|blotch virus|viroid|tobacco|"
    r"tombus|cherry virus|bean endornavirus|coffee ringspot virus|"
    r"Odontoglossum ringspot virus|hibiscus|Obuda pepper virus|"
    r"tobamovirus|"
    r"wheat stripe virus|Lettuce necrotic yellows virus|citrus yellow|"
    r"Maize fine streak virus|dwarf (fiji)?virus|maize stripe virus|"
    r"tomato brown rugose fruit virus|turnip vein-clearing virus|"
    r"chlorotic spot virus)\b",
    re.I,
)

PLANT_ONLY_GENERA = {
    "Emaravirus",  # In Bunyavirales.
    "Idaeovirus",  # Unassigned.
    "Ourmiavirus",  # In Botourmiaviridae.
    "Polemovirus",  # Unassigned.
    "Sobemovirus",  # Unassigned.
    "Tospovirus",  # (Misclassified?) in Peribunyaviridae (e.g., NC_040742.1)
    "Umbravirus",  # In Tombusviridae.
    "Varicosavirus",  # In Rhabdoviridae.
}

PLANT_ONLY_FAMILIES = {
    "Aspiviridae",  # (Formerly Ophioviridae) in Serpentovirales.
    "Avsunviroidae",  # Viroids.
    "Betaflexiviridae",  # In Tymovirales.
    "Bromoviridae",  # Unassigned.
    "Caulimoviridae",  # In Ortervirales.
    "Closteroviridae",  # Unassigned.
    "Geminiviridae",  # Unassigned.
    "Luteoviridae",  # Unassigned.
    "Nanoviridae",  # Unassigned.
    "Partitiviridae",  # Unassigned.
    "Pospiviroidae",  # Unassigned (viroids).
    "Secoviridae",  # In Picornavirales.
    "Tombusviridae",  # Unassigned.
    "Tospoviridae",  # In Bunyavirales.
    "Virgaviridae",  # Unassigned.
}


PLANT_ONLY_ORDERS = {"Tymovirales"}


class _UnknownNameError(Exception):
    "An unknown name was used."


class LineageFetcher:
    """
    Provide access to the NCBI taxonomy database so we can retrieve the lineage
    of title sequences hit by BLAST.
    """

    def __init__(self, db=None, cursor=None):
        self._db = db or getDatabaseConnection()
        self._cursor = cursor or self._db.cursor()
        self._cache = {}

    def lineage(self, title: str) -> list[tuple[int, str]]:
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
        gi = int(title.split("|")[1])
        query = "SELECT taxID FROM gi_taxid WHERE gi = %d" % gi

        try:
            while True:
                self._cursor.execute(query)
                taxID = self._cursor.fetchone()[0]
                if taxID == 1:
                    break
                query = "SELECT name FROM names WHERE taxId = %s" % taxID
                self._cursor.execute(query)
                scientificName = self._cursor.fetchone()[0]
                lineage.append((taxID, scientificName))
                # Move up to the parent.
                query = "SELECT parent_taxID FROM nodes WHERE taxID = %s" % taxID
        except TypeError:
            lineage = []

        self._cache[title] = lineage
        return lineage

    def close(self) -> None:
        """
        Close the database connection and render self invalid. Any subsequent
        re-use of self will raise an error.
        """
        self._cursor.close()
        self._db.close()
        self._cursor = self._db = self._cache = None


class Taxonomy:
    """
    Provide access to the NCBI taxonomy database so we can retrieve the
    taxonomy lineage corresponding to an accession number.

    @param dbFilenameOrConnection: Either a C{str} database filename or an
        open sqlite3 database. The database must contain taxonomy information
        with tables and columns named as in the building scripts used in
        https://github.com/acorg/ncbi-taxonomy-database
    """

    CACHE_SIZE = 10e6

    def __init__(self, dbFilenameOrConnection):
        if isinstance(dbFilenameOrConnection, str):
            self._db = sqliteConnect(dbFilenameOrConnection)
            self._closeConnection = True
        else:
            self._db = dbFilenameOrConnection
            self._closeConnection = False
        self._lineageCache = LRUCache(maxsize=self.CACHE_SIZE)
        self._hostsCache = LRUCache(maxsize=self.CACHE_SIZE)
        self._plantVirusCache = LRUCache(maxsize=self.CACHE_SIZE)
        self._fungusVirusCache = LRUCache(maxsize=self.CACHE_SIZE)

    def lineageFromTaxid(self, taxid, skipFunc=None, stopFunc=None):
        """
        Get lineage information from the taxonomy database for a taxonomy id.

        @param taxid: An C{int} taxonomy id.
        @param skipFunc: A function that takes a C{LineageElement} instance
            and returns C{True} if this lineage item should be excluded from
            the returned C{tuple}.
        @param stopFunc: A function that takes a C{LineageElement} instance
            and returns C{True} if the lineage processing should be stopped.
            The final lineage element is still added to the returned value,
            unless a C{skipFunc} is passed that returns C{True} for it.
        @raise ValueError: If the taxonomy id cannot be found in the names
            table or if any taxonomy id cannot be found in the nodes table
            as we move up the hierarchy.
        @raise AttributeError: If C{self.close} has been called.
        @return: A C{tuple} of the taxonomic categories of the title. Each
            tuple element is a C{LineageElement}, giving a taxonomy id a
            (scientific) name, and the rank (species, genus, etc). The first
            element in the list corresponds to the passed C{taxid} and each
            successive element is the parent of the preceeding one. The top
            level of the taxonomy hierarchy (with taxid = 1) is not included
            in the returned list. Note that the returned C{tuple} may be empty,
            depending on C{skipFunc} and C{stopFunc}.
        """
        cursor = self._db.cursor()
        execute = cursor.execute
        fetchone = cursor.fetchone

        lineage = []

        while taxid != 1:
            execute("SELECT name FROM names WHERE taxid = ?", (taxid,))
            result = fetchone()
            if result is None:
                raise ValueError(
                    "Could not find taxonomy id %r in names table" % (taxid,)
                )
            name = result[0]
            execute("SELECT parent_taxid, rank FROM nodes WHERE taxid = ?", (taxid,))
            result = fetchone()
            if result is None:
                raise ValueError(
                    "Could not find taxonomy id %r in nodes table" % (taxid,)
                )

            parentTaxid, rank = result
            element = LineageElement(taxid, name, rank)

            if not (skipFunc and skipFunc(element)):
                lineage.append(element)

            if stopFunc and stopFunc(element):
                break

            taxid = parentTaxid

        return tuple(lineage)

    @cachedmethod(attrgetter("_lineageCache"))
    def lineage(self, id_, skipFunc=None, stopFunc=None):
        """
        Get lineage information from the taxonomy database for an
        accession number, name, or taxonomy id.

        @param id_: A C{str} accession number, C{str} name (e.g.,
            'Hymenoptera'), or an C{int} taxonomy id. If a C{str} accession
            number, C{id_} must include the version (e.g., the ".1" part of
            "DQ011818.1") if the taxonomy database in use uses it (as is the
            case with the database built by
            https://github.com/acorg/ncbi-taxonomy-database)
        @param skipFunc: A function that takes a C{LineageElement} instance
            and returns C{True} if this lineage item should be excluded from
            the returned C{tuple}.
        @param stopFunc: A function that takes a C{LineageElement} instance
            and returns C{True} if the lineage processing should be stopped.
            The final lineage element is still added to the returned value,
            unless a C{skipFunc} is passed that returns C{True} for it.
        @raise ValueError: If a taxonomy id (or one of its higher-level
            taxonomy ids) cannot be found.
        @raise AttributeError: If C{self.close} has been called.
        @return: A C{tuple} of the taxonomic categories of the title. Each
            tuple element is a C{LineageElement}, giving a taxonomy id a
            (scientific) name, and the rank (species, genus, etc). The first
            element in the tuple corresponds to the passed C{id_} and each
            successive element is the parent of the preceeding one. The top
            level of the taxonomy hierarchy (with taxid = 1) is not included
            in the returned tuple.
        """
        if isinstance(id_, int):
            return self.lineageFromTaxid(id_, skipFunc=skipFunc, stopFunc=stopFunc)

        cursor = self._db.cursor()

        cursor.execute("SELECT taxid FROM accession_taxid WHERE accession = ?", (id_,))
        row = cursor.fetchone()
        if row:
            return self.lineageFromTaxid(
                int(row[0]), skipFunc=skipFunc, stopFunc=stopFunc
            )

        try:
            taxid = self.taxIdFromName(id_)
        except _UnknownNameError:
            raise ValueError(
                "Could not find %r in accession_taxid or names tables" % id_
            )
        else:
            return self.lineageFromTaxid(taxid, skipFunc=skipFunc, stopFunc=stopFunc)

    @cachedmethod(attrgetter("_hostsCache"))
    def hostsFromTaxid(self, taxid):
        """
        Get host information from the taxonomy database for a taxonomy id.

        @param taxid: An C{int} taxonomy id.
        @raise AttributeError: If C{self.close} has been called.
        @return: A C{set} of C{str} host names. If no host information is found
            for C{taxid}, return C{None}.
        """
        cursor = self._db.cursor()
        cursor.execute("SELECT hosts FROM hosts WHERE taxid = ?", (taxid,))
        row = cursor.fetchone()
        if row:
            return set(row[0].split(","))

    def taxIdFromName(self, name):
        """
        Get a taxonomy id from a name.

        @param name: A C{str} taxonomy name. This must be a unique name, as
            found in the names.dmp file of the NCBI taxonomy names.dmp file.
        @raise _UnknownNameError: If a taxonomy id cannot be found for C{name}.
        @return: An C{int} taxonomy id.
        """
        cursor = self._db.cursor()
        cursor.execute("SELECT taxid FROM names WHERE name = ?", (name,))
        row = cursor.fetchone()
        if row:
            return int(row[0])

        raise _UnknownNameError("Name %r not found in names tables" % name)

    def hosts(self, id_):
        """
        Get host information from the taxonomy database for an accession number
        or taxonomy id.

        @param id_: A C{str} accession number or C{int} taxonomy id.
            If a C{str}, C{id_} must include the version (e.g., the
            ".1" part of "DQ011818.1") if the taxonomy database in use uses
            it (as is the case with the database built by
            https://github.com/acorg/ncbi-taxonomy-database
        @raise AttributeError: If C{self.close} has been called.
        @return: A C{set} of C{str} host names, or C{None} if no host
            information.
        """
        if isinstance(id_, int):
            return self.hostsFromTaxid(id_)

        cursor = self._db.cursor()

        cursor.execute("SELECT taxid FROM accession_taxid WHERE accession = ?", (id_,))
        row = cursor.fetchone()
        if row:
            return self.hostsFromTaxid(int(row[0]))

        cursor.execute("SELECT taxid FROM names WHERE name = ?", (id_,))
        row = cursor.fetchone()
        if row:
            return self.hostsFromTaxid(int(row[0]))

    @cachedmethod(attrgetter("_fungusVirusCache"))
    def isFungusOnlyVirus(self, lineage, title=None) -> bool:
        """
        Determine whether a lineage corresponds to a fungus-only virus (i.e.,
        a virus that only infects fungi hosts).

        @param lineage: An iterable of C{tuple}s of taxonomy id, scientific
            name, and rank as returned by C{Taxonomy.lineage}.
        @param title: If not C{None}, a C{str} title of the virus.
        @return: C{True} if the lineage corresponds to a fungus-only virus,
            C{False} otherwise.
        """
        for taxid, name, rank in lineage:
            if (
                rank == "family"
                and name in FUNGUS_ONLY_FAMILIES
                or rank == "genus"
                and name in FUNGUS_ONLY_GENERA
            ):
                return True

        if title and FUNGUS_ONLY_VIRUS_REGEX.search(title):
            return True

        # Do host taxonomy database lookups. Try to look up host
        # information at all taxonomy levels. This is (of course) slower,
        # but it (likely, as with plant viruses) results in more
        # fungus-only viruses being identified.
        for taxid, _, _ in lineage:
            if self.hosts(taxid) == {"fungi"}:
                return True

        return False

    @cachedmethod(attrgetter("_plantVirusCache"))
    def isPlantOnlyVirus(self, lineage, title=None) -> bool:
        """
        Determine whether a lineage corresponds to a plant-only virus (i.e.,
        a virus that only infects plant hosts).

        @param lineage: An iterable of C{tuple}s of taxonomy id, scientific
            name, and rank as returned by C{Taxonomy.lineage}.
        @param title: If not C{None}, a C{str} title of the virus.
        @return: C{True} if the lineage corresponds to a plant-only virus,
            C{False} otherwise.
        """
        for taxid, name, rank in lineage:
            if (
                rank == "order"
                and name in PLANT_ONLY_ORDERS
                or rank == "family"
                and name in PLANT_ONLY_FAMILIES
                or rank == "genus"
                and name in PLANT_ONLY_GENERA
            ):
                return True

        if title and PLANT_ONLY_VIRUS_REGEX.search(title):
            return True

        # Do host taxonomy database lookups. Try to look up host information at
        # all taxonomy levels. This is (of course) slower, but it does result
        # in more plant-only viruses being identified, e.g., with rank
        # 'unclassified' and name 'Partitiviridae'.
        for taxid, _, _ in lineage:
            if self.hosts(taxid) == {"plants"}:
                return True

        return False

    @staticmethod
    def subsetLineageByRanks(lineage, func):
        """
        Extract certain ranks from a lineage.

        @param lineage: An iterable of C{tuple}s of taxonomy id, scientific
            name, and rank as returned by C{Taxonomy.lineage}.
        @param func: A function that accepts a C{str} rank (which may be '-')
            and returns C{True} or C{False} according to whether the lineage
            element should be retained.
        @return: A filter object that yields C{tuple}s of taxonomy id,
            scientific name, and rank, where the elements are those that
            C{func} returns C{True}. The elements in the returned value will
            be in the order in which they are present in C{lineage}.
        """
        return filter(lambda thisLineage: func(thisLineage.rank), lineage)

    def close(self) -> None:
        """
        Close the database connection (if we opened it).
        """
        if self._closeConnection:
            self._db.close()
        # Set self._db to None so self.lineage will raise an AttributeError
        # exception if it is called again.
        self._db = None

    def __enter__(self) -> Taxonomy:
        return self

    def __exit__(self, excType, excValue, traceback) -> None:
        self.close()


def isRetrovirus(lineage: Iterable[LineageElement]) -> bool:
    """
    Determine whether a lineage corresponds to a retrovirus.

    @param lineage: An iterable of C{LineageElement} instances.
    @return: C{True} if the lineage corresponds to a retrovirus, C{False}
        otherwise.
    """
    for element in lineage:
        if element.rank == "family" and element.name == "Retroviridae":
            return True

    return False


def isRNAVirus(lineage):
    """
    Determine whether a lineage corresponds to an RNA virus.

    @param lineage: An iterable of C{LineageElement} instances.
    @return: C{True} if the lineage corresponds to an RNA virus, C{False}
        otherwise.
    """
    for element in lineage:
        if (element.name, element.rank) in RNA_VIRUS_LINEAGE_ELEMENTS:
            return True

    return False


def isDNAVirus(lineage: Iterable[LineageElement]) -> bool:
    """
    Determine whether a lineage corresponds to an DNA virus.

    @param lineage: An iterable of C{LineageElement} instances.
    @return: C{True} if the lineage corresponds to a DNA virus, C{False}
        otherwise.
    """
    return not isRNAVirus(lineage)


def isAllowedTaxonomicRank(allowedTaxonomicRanks, lineage):
    """
    Determine whether a lineage matches a set of desired ranks.

    @param allowedTaxonomicRanks: If not C{None}, a set of (case insensitive)
        acceptable (name, rank) C{str} tuples. E.g.,
            set((
                ("nidovirales", "order"),
                ("retroviridae", "family"),
            ))
    @param lineage: An iterable of C{LineageElement} instances.
    @return: C{True} if the lineage matches at least one of the name:rank
        pairs. Else C{False}.
    """
    lcAllowedTaxonomicRanks = set()
    for name, rank in allowedTaxonomicRanks:
        lcAllowedTaxonomicRanks.add((name.lower(), rank.lower()))

    for element in lineage:
        if (element.name.lower(), element.rank.lower()) in lcAllowedTaxonomicRanks:
            return True

    return False


def formatLineage(
    lineage, namesOnly: bool = False, separator: str = ", ", prefix: str = ""
) -> str:
    """
    Format a lineage for printing.

    @param lineage: An iterable of C{LineageElement} instances.
    @param namesOnly: If C{True} only print taxonomic names.
    @param separator: A C{str} separator to put between fields. If C{None},
        return a space-padded aligned columns.
    @param prefix: A C{str} to put at the start of each line.
    @return: A formatted C{str} for printing.
    """
    if namesOnly:
        # The separator is guaranteed to be set by our caller.
        return prefix + separator.join(element.name for element in lineage)

    if separator is not None:
        return "\n".join(
            "%s%s%s%s%s%s" % (prefix, rank, separator, name, separator, taxid)
            for (taxid, name, rank) in lineage
        )

    # This is a bit slow, walking through the lineage list 3 times.
    taxidWidth = max(len(str(element.taxid)) for element in lineage)
    nameWidth = max(len(element.name) for element in lineage)
    rankWidth = max(len(element.rank) for element in lineage)

    return "\n".join(
        "%s%-*s %-*s %*d"
        % (prefix, rankWidth, rank, nameWidth, name, taxidWidth, taxid)
        for (taxid, name, rank) in lineage
    )


def lineageTaxonomyLinks(lineage: Iterable[LineageElement]) -> list[str]:
    """
    Get HTML links for a lineage.

    @param lineage: An iterable of C{LineageElement} instances.
    @return: A C{list} of HTML C{str} links.
    """
    names = [element.name for element in lineage]
    names[0] = "taxon"

    taxids = [element.taxid for element in lineage]

    return [
        '<a href="%s%s">%s</a>'
        % ("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=", taxid, name)
        for taxid, name in list(zip(taxids, names))[:-1]
    ]


class _HierarchyNode:
    def __init__(self, name: str):
        self.name = name
        self.count = 0
        self.nodes: dict[str, _HierarchyNode] = {}
        self.tips: list[tuple[str, str]] = []

    def addChild(self, name) -> _HierarchyNode:
        self.count += 1
        try:
            node = self.nodes[name]
        except KeyError:
            node = self.nodes[name] = _HierarchyNode(name)

        return node

    def addTip(self, name: str, genomeAccession: str) -> None:
        self.count += 1
        self.tips.append((name, genomeAccession))

    def toDict(self) -> dict:
        nodes = [node.toDict() for node in self.nodes.values()]

        if self.tips:
            nodes.extend(
                [
                    {
                        "text": name,
                        "href": "#pathogen-" + accession,
                    }
                    for (name, accession) in sorted(self.tips)
                ]
            )

        return {
            # 'count': self.count,
            "nodes": nodes,
            # 'name': self.name,
            # 'tips': self.tips,
            "text": "%s (%d)" % (self.name, self.count),
        }


class Hierarchy:
    """
    Collect information about a lineages in a taxonomy hierarchy.
    """

    def __init__(self, rootName: str = "Viruses"):
        self._root = _HierarchyNode(rootName)

    def add(self, lineage: Iterable[LineageElement], genomeAccession: str) -> None:
        """
        Add a new lineage.

        @param lineage: An iterable of C{LineageElement} instances.
        @param genomeAccession: A C{str} pathogen accession number.
        """
        names = [element.name for element in lineage]

        n = len(names) - 2
        node = self._root

        while n:
            node = node.addChild(names[n])
            n -= 1

        node.addTip(names[0], genomeAccession)

    def toJSON(self):
        return dumps([self._root.toDict()])


def addTaxonomyDatabaseCommandLineOptions(parser: argparse.ArgumentParser) -> None:
    """
    Add standard taxonomy database command-line options to an argparse parser.

    @param parser: An C{argparse.ArgumentParser} instance.
    """
    parser.add_argument(
        "--" + TAXONOMY_DATABASE_COMMAND_LINE_OPTION,
        metavar="DATABASE-FILE",
        help=(
            "The file holding an sqlite3 taxonomy database. See "
            "https://github.com/acorg/ncbi-taxonomy-database for how to "
            "build one. If not specified, the value in the %r environment "
            "variable (if any) will be used." % TAXONOMY_DATABASE_ENV_VAR
        ),
    )


def parseTaxonomyDatabaseCommandLineOptions(args, parser) -> Optional[Taxonomy]:
    """
    Examine parsed command-line options and return a Taxonomy instance. Exits
    if no taxonomy database is given or named in the TAXONOMY_DATABASE_ENV_VAR
    environment variable.

    @param args: An argparse namespace, as returned by the argparse
        C{parse_args} function.
    @param parser: An C{argparse.ArgumentParser} instance.
    @return: A C{Taxonomy} instance, or C{None} if no database filename is
        given or set via the environment variable named in
        TAXONOMY_DATABASE_ENV_VAR.
    """
    filename = args.taxonomyDatabase or environ.get(TAXONOMY_DATABASE_ENV_VAR)

    if filename:
        return Taxonomy(filename)
    else:
        parser.error(
            "A taxonomy database file must be given, either via --%s on "
            "the command line or via the %r environment variable."
            % (TAXONOMY_DATABASE_COMMAND_LINE_OPTION, TAXONOMY_DATABASE_ENV_VAR)
        )

        return None

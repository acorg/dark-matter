import matplotlib.pyplot as pyplot
from collections import defaultdict
from dark import mysql


def _getLineageInfoPerTaxID(taxID, cursor):
    """
    For a given taxID, extract information about the lineage from
    a mySQL database.

    @param taxID: a value assigned by NCBI taxonomy which stands
        for taxonomic information.
    @param cursor: A database cursor, with C{execute} and C{fetchone} methods.
    @return: a C{list} containing information about the taxonomic lineage
    """
    result = []
    while taxID != 1:
        questionToNodes = ('SELECT rank, parent_taxID from nodes '
                           'where taxID = %s' % taxID)
        cursor.execute(questionToNodes)
        rank, parentTaxID = cursor.fetchone()[0:2]
        questionToNames = 'SELECT name from names where taxId = %s' % taxID
        cursor.execute(questionToNames)
        scientificName = cursor.fetchone()[0]
        result.append({
            'taxID': taxID,
            'parentTaxID': parentTaxID,
            'rank': rank,
            'scientificName': scientificName,
        })
        taxID = parentTaxID
    return result


def getLineageInfo(blastHits, db=None):
    """
    For each taxID present in blastHits, walk through the taxonomic tree,
    starting from that taxID, and for each node, store the information
    about the node (rank, taxID, parent_taxID, scientificName).

    @param blastHits: A L{dark.blast.BlastHits} instance.
    @param db: A database connection, with C{cursor} and C{close} methods.
    @return: A C{dict} where each key is a taxID and the value is the
        taxonomic information
    """
    if db is None:
        openedDb = True
        db = mysql.getDatabaseConnection()
    else:
        openedDb = False
    cursor = db.cursor()
    taxIDLookUpDict = {}
    for title in blastHits.titles:
        taxID = blastHits.titles[title]['taxID']
        if taxID not in taxIDLookUpDict:
            taxIDLookUpDict[taxID] = _getLineageInfoPerTaxID(taxID, cursor)
    cursor.close()
    if openedDb:
        db.close()
    return taxIDLookUpDict


def taxIDsPerTaxonomicRank(taxIDLookUpDict, taxonomicRank):
    """
    Collects the number of taxIDs that hit at a specified
    taxonomic rank.

    @param taxIDLookUpDict: A C{dict} which contains information
        about the lineage for each taxID. Result of calling
        getLineageInfo().
    @param taxonomicRank: A C{str} specifying the taxonomic rank
        for which taxonomic infromation should be provided. Must
        be one of superkingdom, kingdom, phylum, subphylum,
        superclass, class, superorder, order, suborder, infraorder,
        parvoorder, superfamily, family, subfamily, genus, species.
    @return: A C{dict} where each key is the scientific Name at
        the taxonomic rank specified and the value is the list
        of taxIDs that correspond to that rank.
    """
    taxIDsPerRank = defaultdict(set)
    for daughterTaxID in taxIDLookUpDict:
        for rankItem in taxIDLookUpDict[daughterTaxID]:
            if rankItem['rank'] == taxonomicRank:
                scientificName = rankItem['scientificName']
                taxIDsPerRank[scientificName].add(daughterTaxID)
    return taxIDsPerRank


def readsPerTaxonomicRank(taxIDLookUpDict, blastHits, taxonomicRank):
    """
    Collects the number of reads that hit at a specified
    taxonomic rank.

    @param taxIDLookUpDict: A C{dict} which contains information
        about the lineage for each taxID. Result of calling
        getLineageInfo().
    @param blastHits: A L{dark.blast.BlastHits} instance.
    @param taxonomicRank: A C{str} specifying the taxonomic rank
        for which taxonomic infromation should be provided. Must
        be one of superkingdom, kingdom, phylum, subphylum,
        superclass, class, superorder, order, suborder, infraorder,
        parvoorder, superfamily, family, subfamily, genus, species.
    @return: A C{dict} where each key is the scientific Name at
        the taxonomic rank specified and the value is the list
        of reads that hit at that rank.
    """
    result = defaultdict(set)
    taxIDsPerRank = taxIDsPerTaxonomicRank(taxIDLookUpDict, taxonomicRank)
    for title in blastHits.titles:
        taxID = blastHits.titles[title]['taxID']
        for rank, taxIDs in taxIDsPerRank.iteritems():
            if taxID in taxIDs:
                for item in blastHits.titles[title]['plotInfo']['items']:
                    if item['readNum'] not in result[rank]:
                        result[rank].add(item['readNum'])
    return result


def subjectsPerTaxonomicRank(taxIDLookUpDict, blastHits, taxonomicRank):
    """
    Collects the number of subjects that are hit at a specified
    taxonomic rank.

    @param taxIDLookUpDict: A C{dict} which contains information
        about the lineage for each taxID. Result of calling
        getLineageInfo().
    @param blastHits: A L{dark.blast.BlastHits} instance.
    @param taxonomicRank: A C{str} specifying the taxonomic rank
        for which taxonomic infromation should be provided. Must
        be one of superkingdom, kingdom, phylum, subphylum,
        superclass, class, superorder, order, suborder, infraorder,
        parvoorder, superfamily, family, subfamily, genus, species.
    @return: A C{dict} where each key is the scientific Name at
        the taxonomic rank specified and the value is the list
        of subjects that correspond to that rank.
    """
    result = defaultdict(set)
    taxIDsPerRank = taxIDsPerTaxonomicRank(taxIDLookUpDict, taxonomicRank)
    for title in blastHits.titles:
        taxID = blastHits.titles[title]['taxID']
        for rank, taxIDs in taxIDsPerRank.iteritems():
            if taxID in taxIDs:
                result[rank].add(title)
    return result


def pieChartPerSpecificRank(blastHits, taxonomicRank, what):
    """
    Makes a pie chart with the information specified.

    @param blastHits: A L{dark.blast.BlastHits} instance.
    @param taxonomicRank: A C{str} specifying the taxonomic rank
        for which taxonomic infromation should be provided. Must
        be one of superkingdom, kingdom, phylum, subphylum,
        superclass, class, superorder, order, suborder, infraorder,
        parvoorder, superfamily, family, subfamily, genus, species.
    @param what: A C{str} specifying what information should be
        plotted. Must be one of 'subjects', 'reads', 'taxIDs'.
    """
    taxIDLookUpDict = getLineageInfo(blastHits)
    assert what in ('subjects', 'reads', 'taxIDs'), (
        'what must be one of "subjects", "reads", "taxIDs"')

    assert taxonomicRank in ('superkingdom', 'kingdom', 'phylum', 'subphylum',
                             'superclass', 'class', 'superorder', 'order',
                             'suborder', 'infraorder', 'parvoorder',
                             'superfamily', 'family', 'subfamily', 'genus',
                             'species'), (
        'taxonomicRank must be one of "superkingdom", "kingdom", "phylum", '
        '"subphylum", "superclass", "class", "superorder", "order", '
        '"suborder", "infraorder", "parvoorder", "superfamily", "family", '
        '"subfamily", "genus", "species"')

    if what == 'subjects':
        result = subjectsPerTaxonomicRank(taxIDLookUpDict, blastHits,
                                          taxonomicRank)
    elif what == 'reads':
        result = readsPerTaxonomicRank(taxIDLookUpDict, blastHits,
                                       taxonomicRank)
    elif what == 'taxIDs':
        result = taxIDsPerTaxonomicRank(taxIDLookUpDict, taxonomicRank)

    toPlot = []
    legend = []
    for rank, read in result.iteritems():
        toPlot.append(len(read))
        title = '%s, %s' % (len(read), rank)
        legend.append(title)
    plotTitle = 'Number of %s per %s' % (what, taxonomicRank)
    pyplot.axis('equal')
    pyplot.pie(toPlot)
    pyplot.title(plotTitle, fontsize=20)
    pyplot.legend(legend, bbox_to_anchor=(2, 1.05))
    pyplot.show()


class LineageFetcher(object):
    """
    Provide access to the NCBI taxonomy database so we can retrieve the lineage
    of title sequences hit by BLAST.
    """
    def __init__(self):
        self._db = mysql.getDatabaseConnection()
        self._cursor = self._db.cursor()
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
        query = 'SELECT taxID from gi_taxid_nucl where gi = %d' % gi
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

import matplotlib.pyplot as pyplot
from collections import defaultdict
from dark import mysql


def _getLineageInfoPerTaxID(taxID, db):
    """
    For a given taxID, extract information about the lineage from
    a mySQL database.

    @param taxID: a value assigned by NCBI taxonomy which stands
        for taxonomic information.
    @return: a C{list} containing information about the taxonomic lineage
    """
    result = []
    # TODO: the arguments for connecting to the database should be
    # changed once the setup is finalized.
    cursor = db.cursor()
    while taxID != 1:
        questionToNodes = ('SELECT rank, parent_taxID from nodes '
                           'where taxID = %d' % taxID)
        cursor.execute(questionToNodes)
        rank, parentTaxID = cursor.fetchone()[0:2]
        questionToNames = 'SELECT name from names where taxId = %d' % taxID
        cursor.execute(questionToNames)
        scientificName = cursor.fetchone()[0]
        result.append({
            'taxID': taxID,
            'parentTaxID': parentTaxID,
            'rank': rank,
            'scientificName': scientificName,
        })
        taxID = parentTaxID
    cursor.close()
    return result


def getLineageInfo(blastHits, db=None):
    """
    For each taxID present in blastHits, walk through the taxonomic tree,
    starting from that taxID, and for each node, store the information
    about the node (rank, taxID, parent_taxID, scientificName).

    @param blastHits: A L{dark.blast.BlastHits} instance.
    @return: A C{dict} where each key is a taxID and the value is the
        taxonomic information
    """
    db = db or mysql.getDatabaseConnection()
    taxIDLookUpDict = {}
    for title in blastHits.titles:
        taxID = blastHits.titles[title]['taxID']
        if taxID not in taxIDLookUpDict:
            taxIDLookUpDict[taxID] = _getLineageInfoPerTaxID(taxID, db)
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

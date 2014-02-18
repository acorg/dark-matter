import MySQLdb
import matplotlib.pyplot as pyplot
from collections import defaultdict


def _getLineageInfoPerTaxID(taxID):
    """
    For a given taxID, extract information about the lineage from
    a mySQL database.

    @param taxID: a value assigned by NCBI taxonomy which stands
        for taxonomic information.
    """
    result = []
    # TODO: the arguments for connecting to the database should be
    # changed once the setup is finalized.
    db = MySQLdb.connect(host='localhost', user='root',
                         passwd='rootpassword', db='ncbi_taxonomy')
    cursor = db.cursor()
    n = 1
    while taxID != 1:
        questionToNodes = ('SELECT rank, parent_taxID from nodes '
                           'where taxID = %d' % taxID)
        cursor.execute(questionToNodes)
        rankParent = cursor.fetchone()
        thisResult = {
            'index': n,
            'taxID': taxID,
            'parentTaxID': rankParent[1],
            'rank': rankParent[0],
            'scientificName': None,
        }
        questionToNames = 'SELECT name from names where taxId = %d' % taxID
        cursor.execute(questionToNames)
        names = cursor.fetchone()
        thisResult['scientificName'] = names[0]
        result.append(thisResult)
        taxID = thisResult['parentTaxID']
        n += 1
    cursor.close()
    db.close()
    return result


def getLineageInfo(blastHits):
    """
    For each taxID present in blastHits, walk through the taxonomic tree,
    starting from that taxID, and for each node, store the information
    about the node (rank, taxID, parent_taxID, scientificName).

    @param blastHits: A L{dark.blast.BlastHits} instance.
    """
    taxIDLookUpDict = {}
    for title in blastHits.titles:
        if blastHits.titles[title]['taxID'] not in taxIDLookUpDict:
            taxID = blastHits.titles[title]['taxID']
            lineage = _getLineageInfoPerTaxID(taxID)
            taxIDLookUpDict[blastHits.titles[title]['taxID']] = lineage
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
    """
    taxIDsPerRank = defaultdict(list)
    for daughterTaxID in taxIDLookUpDict:
        for rankItem in taxIDLookUpDict[daughterTaxID]:
            if rankItem['rank'] == taxonomicRank:
                scientificName = rankItem['scientificName']
                taxIDsPerRank[scientificName].append(daughterTaxID)
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
    """
    result = defaultdict(list)
    taxIDsPerRank = taxIDsPerTaxonomicRank(taxIDLookUpDict, taxonomicRank)
    for title in blastHits.titles:
        taxID = blastHits.titles[title]['taxID']
        for rank, taxIDs in taxIDsPerRank.iteritems():
            if taxID in taxIDs:
                for item in blastHits.titles[title]['plotInfo']['items']:
                    if item['readNum'] not in result[rank]:
                        result[rank].append(item['readNum'])
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
    """
    result = defaultdict(list)
    taxIDsPerRank = taxIDsPerTaxonomicRank(taxIDLookUpDict, taxonomicRank)
    for title in blastHits.titles:
        taxID = blastHits.titles[title]['taxID']
        for rank, taxIDs in taxIDsPerRank.iteritems():
            if taxID in taxIDs:
                result[rank].append(title)
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
    print legend
    plotTitle = 'Number of %s per %s' % (what, taxonomicRank)
    pyplot.axis('equal')
    pyplot.pie(toPlot)
    pyplot.title(plotTitle, fontsize=20)
    pyplot.legend(legend, bbox_to_anchor=(2, 1.05))
    pyplot.show()

import matplotlib.pylab as plt
import numpy as np
import re
from scipy.cluster.vq import kmeans, vq
from Bio import SeqIO
import sys

from dark import blast
from dark.dimension import dimensionalIterator


GULLREGEX = re.compile('gull|kittiwake', re.I)
DUCKREGEX = re.compile('duck|mallard|teal|garganey|gadwall|shoveler|'
                       'pochard|wigeon|merganser|nene|scaup|screamer|'
                       'goose|swan|pintail|smew|scoter|eider|bufflehead|'
                       'goldeneye', re.I)

NEITHER = 0
GULL = 1
DUCK = 2


# General utility functions

def _getBird(title):
    """
    Finds out whether a title of a sequence belongs to a bird or a gull.

    @param title: The title of a blastHit.
    """
    gull = GULLREGEX.search(title)
    duck = DUCKREGEX.search(title)

    if gull:
        return GULL
    elif duck:
        return DUCK
    else:
        return NEITHER


def computePercentId(alignment):
    """
    Calculates the percent sequence identity of an alignment.

    @param alignment: An instance of C{Bio.Blast.Record.Blast.Alignment}
    """
    query = alignment.hsps[0].query
    sbjct = alignment.hsps[0].sbjct
    length = len(query)

    identical = 0
    for queryBase, subjectBase in zip(query, sbjct):
        identical += int(queryBase == subjectBase)

    identity = identical / float(length) * 100
    return identity


def makeListOfHitTitles(blastName):
    """
    Makes a list of titles that are hit.

    @param blastName: File with blast output
    """
    blastRecords = blast.BlastRecords(blastName)
    interesting = blastRecords.filterHits(withBitBetterThan=50)
    titlesList = interesting.titles.keys()

    return titlesList


# functions for working with distance graphs

def getPositionInMatrix(matrix, title):
    """
    Given a matrix with titles, return possible coordinates of
    the matrix.

    @matrix: A distance matrix with borders, as returned by
        makeDistanceMatrix().
    @title: A C{str} of a sequence title.
    """
    y = range(len(matrix[0]))
    x = range(len(matrix))

    for xs in x:
        for ys in y:
            if matrix[xs][ys] == title:
                return xs, ys
    print >>sys.stderr, '%s could not be found.' % title
    return False, False


def distancePlot(record, readsAx=False, distance='bit'):
    """
    Produces a rectangular panel of graphs that each show sorted distances for
    a read. Read hits against a certain strain (see find, below) are
    highlighted.

    @param record: A C{Bio.Blast.Record.Blast} instance.
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param distance: The measure of distance read out from the blastFile,
        either 'bit' of 'percentId'.

    @return: Returns the largest distance, and the number of distances that
        were plotted.
    """
    alignments = record.alignments
    title = str(record.query)

    fig = plt.figure(figsize=(10, 10))
    ax = readsAx or fig.add_subplot(111)

    sortedAlignments = sorted(alignments,
                              key=lambda k: k.hsps[0].bits,
                              reverse=True)

    distances = []
    titles = []

    for alignment in sortedAlignments:
        if distance == 'bit':
            distances.append(alignment.hsps[0].bits)
        else:
            dist = computePercentId(alignment)
            distances.append(dist)
        titles.append(alignment.title)

    x = np.arange(0, len(distances))

    # plot black line with distance
    ax.plot(x, distances, 'k')

    # plot green and red dots denoting gulls and ducks
    for i, alignment in enumerate(sortedAlignments):
        title = alignment.title
        y = distances[i]
        bird = _getBird(title)
        if bird == 'gull':
            ax.plot([i], [y], 'rx')
        elif bird == 'duck':
            ax.plot([i], [y], 'gx')

    return distances[0], len(distances)


def distancePanel(blastName, matrix, distance='bit'):
    """
    Make a panel of distance plots generated with the distancePlot
    function above.

    @param blastName: File with blast output
    @param matrix: A matrix of strings corresponding to record.queries
        at the position where the plot of a given record should be.
    @param distance: The measure of distance read out from the blastFile,
        either 'bit' or 'percentId'.
    """
    cols = 8
    rows = 53
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    maxDistance = 0
    maxReads = 0
    blastRecords = blast.BlastRecords(blastName)
    records = blastRecords.records()

    for record in records:
        query = record.query
        try:
            row, col = getPositionInMatrix(matrix, query)
        except TypeError:
            # if that record is not present in matrix, leave it out.
            continue
        localMaxDistance, numberOfReads = distancePlot(
            record, readsAx=ax[row][col])
        try:
            title = query.split('(')[1]
            subtype = query.split('(')[2][:-2]
            segment = query.split(' ')[0]
        except IndexError:
            # this is for one title which doesn't fit the usual format
            segment = 'Segment 2'
            subtype = 'H3N8'
            title = query
        ax[row][col].set_title('%s, %d, %s \n %s' % (
                               segment, numberOfReads,
                               subtype, title), fontsize=10)
        if localMaxDistance > maxDistance:
            maxDistance = localMaxDistance
        if numberOfReads > maxReads:
            maxReads = numberOfReads

    coords = dimensionalIterator((rows, cols))

    blastRecords = blast.BlastRecords(blastName)
    records = blastRecords.records()

    # normalize plots
    for i, record in enumerate(records):
        row, col = coords.next()
        a = ax[row][col]
        a.set_ylim([0, maxDistance + 1])
        a.axhline(1, color='k')
        a.axhline(maxDistance + 1, color='k')
        a.set_xlim([0, maxReads + 1])
        a.set_yticks([])
        a.set_xticks([])

    for row, col in coords:
        ax[row][col].axis('off')

    plt.subplots_adjust(hspace=0.4)
    figure.suptitle('X: 0 to %d, Y (%s): 0 to %d' %
                    (maxReads, ('Bit scores' if distance == 'bit'
                                else '%\ id'),
                     maxDistance), fontsize=20)
    figure.set_size_inches(5 * cols, 3 * rows, forward=True)


# Functions for working with distance matrix

def makeDistanceMatrix(blastName, fastaName, titlesList=False,
                       masked=False, toFile=False, addTitles=False,
                       missingValue=0, distance='bit'):
    """
    Takes a blast output file, returns a distance matrix.

    @param blastName: File with blast output
    @param fastaName: A fastafile with the titles that were blasted,
        in the order that it should be in the matrix. Or a list of
        titles that were blasted, in the order that they should be
        in the matrix.
    @param titlesList: If not C{False}, a list of titles, in the order
        that they should be in the matrix.
    @param masked: If C{True}, the matrix returned will have cells with no hits
        masked out.
    @param toFile: If not C{False}, a C{str} file name where thematrix should
        be written to.
    @param addTitles: If C{True} the titles of hits and fasta sequences will be
        added to the matrix.
    @param missingValue: The value used in the matrix if no distance is given.
    @param distance: The measure of distance read out from the blastFile,
        either 'bit' or 'percentId'.

    NOTE: masked = True and missingValue != 0 does not work!
    """
    if not titlesList:
        titlesList = makeListOfHitTitles(blastName)
    titlesDict = {item: index for (index, item) in enumerate(titlesList)}

    if type(fastaName) == list:
        fastaList = fastaName
    else:
        fastaList = list(record.description for record
                         in SeqIO.parse(fastaName, 'fasta'))
    fastaDict = {item: index for (index, item) in enumerate(fastaList)}

    if not masked:
        initMatrix = np.array(
            [[missingValue for _ in range(
                len(titlesList))] for _ in range(len(fastaList))])
    else:
        initMatrix = np.mp.masked_array(
            [[missingValue for _ in range(
                len(titlesList))] for _ in range(len(fastaList))])

    BlastRecords = blast.BlastRecords(blastName)
    records = BlastRecords.records(bitScoreCutoff=50)

    for record in records:
        query = record.query
        # get position of query in queryList
        try:
            queryIndex = fastaDict[query]
        except ValueError:
            # if query not present in fastaDict, continue
            continue
        for alignment in record.alignments:
            title = alignment.title
            if distance == 'bit':
                dist = alignment.hsps[0].bits
            else:
                dist = computePercentId(alignment)
            # get position of title in titlesList
            try:
                subjectIndex = titlesDict[title]
            except ValueError:
                # if query not present in titlesDict, continue
                continue
            # add distance to matrix
            initMatrix[queryIndex][subjectIndex] = dist

    # mask values that are 0 (those that were not hit)
    if masked:
        initMatrix = np.ma.masked_array(initMatrix, mask=(initMatrix < 1))

    if addTitles:
        initMatrix = distanceMatrixWithBorders(initMatrix,
                                               titlesList, fastaList)

    if toFile:
        matrixToFile(toFile, initMatrix)

    return initMatrix, titlesList, fastaList


def distanceMatrixWithBorders(matrix, titlesList, fastaList):
    """
    Add titles to distance matrix.

    @param matrix: Matrix as returned from makeDistanceMatrix.
    @param titlesList: A list of titles, in the order that they
        should be in the matrix
    @param fastaList: List of titles that were blasted,
        in the order that they should be in the matrix
    """
    for matrixElement, fastaElement in zip(matrix, fastaList):
        matrixElement.insert(0, fastaElement)

    # add an empty element to the beginning of titlesList.
    titlesList.insert(0, '')
    # insert titlesList as the first element of the matrix.
    matrix.insert(0, titlesList)

    return matrix


def matrixToFile(fileName, matrix):
    """
    Writes a matrix to a tab delimited file.

    @param fileName: The name of the file where the distance matrix
        should be written to.
    @param matrix: A distance matrix, as returned from makeDistanceMatrix.
    """
    first = True
    with open(fileName, 'w') as fp:
        for item in matrix:
            if first:
                first = False
                continue
            else:
                fp.write('\n')
            for nr in item:
                fp.write(str(nr) + '\t')


def kMeansCluster(matrix, k):
    """
    Takes a matrix, runs k-means clustering
    Adapted from http://glowingpython.blogspot.co.uk/
    2012/04/k-means-clustering-with-scipy.html

    @param matrix: A matrix (numpy array) for clustering.
    @param k: Number of clusters.
    """
    # computing k-means
    centroids, distortion = kmeans(matrix, 2)

    # assign each sample to a cluster
    index, distortion = vq(matrix, centroids)

    # plotting
    plt.plot(matrix[index == 0, 0], matrix[index == 0, 1], 'bo',
             matrix[index == 1, 0], matrix[index == 1, 1], 'ro')
    plt.plot(centroids[:, 0], centroids[:, 1], 'gs', markersize=8)
    plt.show()


def _annotateTitles(titlesList):
    """
    Takes a list of titles, and returns a list of tuples,
    where each tuple is a title and 0, 1, 2, depending
    on whether the title is a gull or a duck.

    @param titlesList: A list of titles.
    """
    annotatedTitlesList = [(title, _getBird(title)) for title in titlesList]

    return annotatedTitlesList


def distancesBoxPlot(blastName, fastaName, plotTitle, distance='bit'):
    """
    Draws a boxplot of the distances between ducks and gulls.

    @param blastName: file with blast output
    @param fastaName: A fastafile with the titles that was blasted,
        in the order that it should be in the matrix.
    @param plotTitle: a C{str} title of the plot
    """
    matrix, titlesList, fastaList = makeDistanceMatrix(blastName, fastaName,
                                                       missingValue=0,
                                                       distance=distance)

    annotatedTitlesList = _annotateTitles(titlesList)
    annotatedFastaList = _annotateTitles(fastaList)

    gullgullDistances = []
    duckduckDistances = []
    duckgullDistances = []
    gullduckDistances = []

    rows = len(matrix) - 1
    cols = len(matrix[0]) - 1
    for row in range(rows):
        for col in range(cols):
            # skip if no bird was assigned
            if annotatedFastaList[row][1] == 0:
                continue
            if matrix[row][col] != 0:
                if (annotatedFastaList[row][1] == 1 and
                        annotatedTitlesList[col][1] == 1):
                    gullgullDistances.append(matrix[row][col])
                elif (annotatedFastaList[row][1] == 1 and
                      annotatedTitlesList[col][1] == 2):
                    gullduckDistances.append(matrix[row][col])
                elif (annotatedFastaList[row][1] == 2 and
                      annotatedTitlesList[col][1] == 2):
                    duckduckDistances.append(matrix[row][col])
                elif (annotatedFastaList[row][1] == 2 and
                      annotatedTitlesList[col][1] == 1):
                    duckgullDistances.append(matrix[row][col])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    distances = [gullgullDistances, duckduckDistances,
                 duckgullDistances, gullduckDistances]

    ax.boxplot(distances)
    ax.set_xticklabels(['Gull-Gull', 'Duck-Duck', 'Duck-Gull', 'Gull-Duck'])
    ax.set_title(plotTitle)
    ax.set_ylabel('Bit score' if distance == 'bit' else '%\ id')

    return (gullgullDistances, duckduckDistances,
            duckgullDistances, gullduckDistances)

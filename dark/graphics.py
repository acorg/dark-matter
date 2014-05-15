import os
from stat import S_ISDIR
from math import ceil
from collections import defaultdict
from time import ctime, time
from math import log10
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec
from IPython.display import HTML
import numpy as np
from operator import itemgetter

from dark.baseimage import BaseImage
from dark.dimension import dimensionalIterator
from dark.html import AlignmentPanelHTML, NCBISequenceLink, NCBISequenceLinkURL
from dark.intervals import ReadIntervals
from dark.features import ProteinFeatureAdder, NucleotideFeatureAdder
from dark import ncbidb
from dark import orfs


QUERY_COLORS = {
    'A': (1.0, 0.0, 0.0),  # Red.
    'C': (0.0, 0.0, 1.0),  # Blue.
    'G': (0.0, 1.0, 0.0),  # Green.
    'N': (1.0, 0.0, 1.0),  # Purple.
    'T': (1.0, 0.8, 0.0),  # Orange.
    'gap': (0.2, 0.2, 0.2),  # Almost black.
    'match': (0.9, 0.9, 0.9),  # Almost white.
    '*': (0.9, 0.9, 0.9),  # Almost white.
}

DEFAULT_BASE_COLOR = (0.5, 0.5, 0.5)  # Grey

# If we're making a plot that has a log-linear X axis, don't show
# background light grey rectangles for any gap whose (logged) width is less
# than SMALLEST_LOGGED_GAP_TO_DISPLAY.
SMALLEST_LOGGED_GAP_TO_DISPLAY = 20


def report(msg):
    print '%s: %s' % (ctime(time()), msg)


def _sortHTML(hits, by):
    """
    Return an C{IPython.display.HTML} object with the hits sorted by the
    given attribute.

    @param hits: An L{dark.blast.BlastHits} instance.
    @param by: A C{str}, one of 'eMean', 'eMedian', 'readCount', 'title'.
    @return: An HTML instance with sorted titles and information about
        hit read count, length, and e-values.
    """
    out = []
    for i, title in enumerate(hits.sortTitles(by), start=1):
        hitInfo = hits.titles[title]
        link = NCBISequenceLink(title, title)
        out.append(
            '%3d: count=%4d, len=%7d, min(bit)=%20s median(bit)=%20s: %s' %
            (i, hitInfo['readCount'], hitInfo['length'],
             hitInfo['minBitScore'], hitInfo['medianBitScore'], link))
    return HTML('<pre><tt>' + '<br/>'.join(out) + '</tt></pre>')


def summarizeHitsByMeanEValue(hits):
    """
    Sort hit titles by mean e-value.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        mean e-value.
    """
    return _sortHTML(hits, 'eMean')


def summarizeHitsByMedianEValue(hits):
    """
    Sort hit titles by median e-value.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        median e-value.
    """
    return _sortHTML(hits, 'eMedian')


def summarizeHitsByMinEValue(hits):
    """
    Sort hit titles by minimum (i.e., best) e-value.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        minimum e-value.
    """
    return _sortHTML(hits, 'eMedian')


def summarizeHitsByCount(hits):
    """
    Sort hit titles by read count.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        read count.
    """
    return _sortHTML(hits, 'readCount')


def summarizeHitsByLength(hits):
    """
    Sort hit titles by sequence length.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        mean sequence length.
    """
    return _sortHTML(hits, 'length')


def summarizeHitsByMaxBitScore(hits):
    """
    Sort hit titles by sequence length.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        max bitscore.
    """
    return _sortHTML(hits, 'bitScoreMax')


def summarizeHitsByMeanBitScore(hits):
    """
    Sort hit titles by sequence length.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        mean bitScore.
    """
    return _sortHTML(hits, 'bitScoreMean')


def summarizeHitsByMedianBitScore(hits):
    """
    Sort hit titles by sequence length.

    @param hits: An L{dark.blast.BlastHits} instance.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        median bit score.
    """
    return _sortHTML(hits, 'bitScoreMedian')


def alignmentGraph(blastHits, title, addQueryLines=True, showFeatures=True,
                   colorQueryBases=False, createFigure=True, showFigure=True,
                   readsAx=None, imageFile=None, quiet=False, idList=False,
                   xRange='subject', plot='e values', showOrfs=True):
    """
    Align a set of matching reads against a BLAST hit.

    @param blastHits: A L{dark.blast.BlastHits} instance.
    @param title: A C{str} sequence title that was hit by BLAST. We plot the
        reads that hit this title.
    @param addQueryLines: if C{True}, draw query lines in full (these will then
        be partly overdrawn by the HSP match against the subject). These are
        the 'whiskers' that potentially protrude from each side of a query.
    @param showFeatures: if C{True}, look online for features of the subject
        sequence (given by hitId).
    @param colorQueryBases: if C{True}, color each base of a query string. If
        C{True}, then addQueryLines is meaningless since the whole query is
        shown colored.
    @param createFigure: If C{True}, create a figure and give it a title.
    @param showFigure: If C{True}, show the created figure. Set this to
        C{False} if you're creating a panel of figures or just want to save an
        image (with C{imageFile}).
    @param readsAx: If not None, use this as the subplot for displaying reads.
    @param imageFile: If not None, specifies a filename to write the image to.
    @param quiet: If C{True}, don't print progress / timing output.
    @param idList: a dictionary. The keys is a color and the values is a list
        of read identifiers that should be colored in the respective color.
    @param xRange: set to either 'subject' or 'reads' to indicate the range of
        the X axis.
    @param plot: A C{str}, either 'bit scores' or 'e values', to indicate
        which values should be plotted on the graph Y axis.
    @param showOrfs: If C{True}, open reading frames will be displayed.
    """

    start = time()

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    assert plot in ('bit scores', 'e values'), (
        'plot must be either "bit scores", or "e values".')

    params = blastHits.plotParams
    assert params is not None, ('Oops, it looks like you forgot to run '
                                'computePlotInfo.')

    application = blastHits.records.blastParams['application'].lower()

    if application == 'blastx' and showOrfs:
        # We cannot show ORFs when displaying protein plots.
        showOrfs = False

    sequence = ncbidb.getSequence(title, blastHits.records.blastDb)

    if createFigure:
        width = 20
        figure = plt.figure(figsize=(width, 20))

    createdReadsAx = readsAx is None

    if showFeatures:
        if showOrfs:
            gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 12])
            featureAx = plt.subplot(gs[0, 0])
            orfAx = plt.subplot(gs[1, 0])
            orfReversedAx = plt.subplot(gs[2, 0])
            readsAx = readsAx or plt.subplot(gs[3, 0])
        else:
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
            featureAx = plt.subplot(gs[0, 0])
            readsAx = readsAx or plt.subplot(gs[1, 0])
    else:
        if showOrfs:
            gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 12])
            orfAx = plt.subplot(gs[0, 0])
            orfReversedAx = plt.subplot(gs[1, 0])
            readsAx = readsAx or plt.subplot(gs[2, 0])
        else:
            readsAx = readsAx or plt.subplot(111)

    plotInfo = blastHits.titles[title]['plotInfo']

    logLinearXAxis = params['logLinearXAxis']

    if logLinearXAxis:
        offsetAdjuster = plotInfo['offsetAdjuster'].adjustOffset
    else:
        offsetAdjuster = lambda x: x

    items = plotInfo['items']
    if plot == 'bit scores':
        maxYIncludingRandoms = maxY = int(ceil(plotInfo['bitScoreMax']))
        minY = int(plotInfo['bitScoreMin'])
    else:
        maxYIncludingRandoms = int(ceil(plotInfo['maxEIncludingRandoms']))
        maxY = int(ceil(plotInfo['maxE']))
        minY = int(plotInfo['minE'])
    maxX = plotInfo['maxX']
    minX = plotInfo['minX']

    # Add light grey vertical rectangles to show the logarithmic gaps. Add
    # these first so that reads will be plotted on top of them. Only draw
    # gaps that are more than SMALLEST_LOGGED_GAP_TO_DISPLAY pixels wide as
    # we could have millions of tiny gaps for a bacteria and drawing them
    # all will be slow and only serves to make the entire background grey.
    if logLinearXAxis and len(plotInfo['offsetAdjuster'].adjustments()) < 100:
        for (intervalType, interval) in plotInfo['readIntervals'].walk():
            if intervalType == ReadIntervals.EMPTY:
                adjustedStart = offsetAdjuster(interval[0])
                adjustedStop = offsetAdjuster(interval[1])
                width = adjustedStop - adjustedStart
                if width >= SMALLEST_LOGGED_GAP_TO_DISPLAY:
                    readsAx.axvspan(adjustedStart, adjustedStop,
                                    color='#f4f4f4')

    # A function to pull the Y value out of a plotInfo.
    getY = itemgetter('bitScore' if plot == 'bit scores' else 'convertedE')

    if colorQueryBases:
        # Color each query by its bases.
        xScale = 3
        yScale = 2
        baseImage = BaseImage(
            maxX - minX,
            maxYIncludingRandoms - minY + (1 if params['rankValues'] else 0),
            xScale, yScale)
        for item in items:
            hsp = item['hsp']
            y = getY(item) - minY
            # If the product of the subject and query frame values is +ve,
            # then they're either both +ve or both -ve, so we just use the
            # query as is. Otherwise, we need to reverse complement it.
            if item['frame']['subject'] * item['frame']['query'] > 0:
                query = blastHits.fasta[item['readNum']].seq
            else:
                # One of the subject or query has negative sense.
                query = blastHits.fasta[
                    item['readNum']].reverse_complement().seq
            query = query.upper()
            queryStart = hsp['queryStart']
            # There are 3 parts of the query string we need to display. 1)
            # the left part (if any) before the matched part of the
            # subject. 2) the matched part (which can include gaps in the
            # query and/or subject). 3) the right part (if any) after the
            # matched part.
            # For each part, calculate the ranges in which we have to make
            # the comparison between subject and query.

            # NOTE: never use item['origHsp'].gaps to calculate the number
            # of gaps, as this number contains gaps in both subject and
            # query.

            # 1. Left part:
            leftRange = hsp['subjectStart'] - queryStart

            # 2. Match, middle part:
            middleRange = len(item['origHsp'].query)

            # 3. Right part:
            # Using hsp['queryEnd'] - hsp['subjectEnd'] to calculate the length
            # of the right part leads to the part being too long. The number of
            # gaps needs to be subtracted to get the right length.
            origQuery = item['origHsp'].query
            rightRange = hsp['queryEnd']-hsp['subjectEnd']-origQuery.count('-')

            # 1. Left part.
            xOffset = queryStart - minX
            queryOffset = 0
            for queryIndex in xrange(leftRange):
                color = QUERY_COLORS.get(query[queryOffset + queryIndex],
                                         DEFAULT_BASE_COLOR)
                baseImage.set(xOffset + queryIndex, y, color)

            # 2. Match part.
            xOffset = hsp['subjectStart'] - minX
            xIndex = 0
            queryOffset = hsp['subjectStart'] - hsp['queryStart']
            origSubject = item['origHsp'].sbjct
            origQuery = item['origHsp'].query.upper()
            for matchIndex in xrange(middleRange):
                if origSubject[matchIndex] == '-':
                    # A gap in the subject was needed to match the query.
                    # In our graph we keep the subject the same even in the
                    # case where BLAST opened gaps in it, so we compensate
                    # for the gap in the subject by not showing this base
                    # of the query.
                    pass
                else:
                    if origSubject[matchIndex] == origQuery[matchIndex]:
                        # The query matched the subject at this location.
                        # Matching bases are all colored in the same
                        # 'match' color.
                        color = QUERY_COLORS['match']
                    else:
                        if origQuery[matchIndex] == '-':
                            # A gap in the query. All query gaps get the
                            # same 'gap' color.
                            color = QUERY_COLORS['gap']
                        else:
                            # Query doesn't match subject (and is not a gap).
                            color = QUERY_COLORS.get(origQuery[matchIndex],
                                                     DEFAULT_BASE_COLOR)
                    baseImage.set(xOffset + xIndex, y, color)
                    xIndex += 1

            # 3. Right part.
            xOffset = hsp['subjectEnd'] - minX
            backQuery = query[-rightRange:].upper()
            for queryIndex in xrange(rightRange):
                color = QUERY_COLORS.get(backQuery[queryIndex],
                                         DEFAULT_BASE_COLOR)
                baseImage.set(xOffset + queryIndex, y, color)

        readsAx.imshow(baseImage.data, aspect='auto', origin='lower',
                       interpolation='nearest',
                       extent=[minX, maxX, minY, maxYIncludingRandoms])
    else:
        # Add horizontal lines for all the query sequences. These will be the
        # grey 'whiskers' in the plots once we (below) draw the matched part
        # on top of part of them.
        if addQueryLines:
            for item in items:
                y = getY(item)
                hsp = item['hsp']
                line = Line2D([hsp['queryStart'], hsp['queryEnd']], [y, y],
                              color='#aaaaaa')
                readsAx.add_line(line)

        # Add the horizontal BLAST alignment lines.

        # If an idList is given set things up to look up read colors.
        readColor = {}
        if idList:
            for color, reads in idList.iteritems():
                for read in reads:
                    if read in readColor:
                        raise ValueError('Read %s is specified multiple '
                                         'times in idList' % read)
                    else:
                        readColor[read] = color

        for item in items:
            queryId = blastHits.fasta[item['readNum']].id
            y = getY(item)
            hsp = item['hsp']
            line = Line2D([hsp['subjectStart'], hsp['subjectEnd']], [y, y],
                          color=readColor.get(queryId, 'blue'))
            readsAx.add_line(line)

    if showOrfs:
        orfs.addORFs(orfAx, sequence.seq, minX, maxX, offsetAdjuster)
        orfs.addReversedORFs(orfReversedAx,
                             sequence.reverse_complement().seq,
                             minX, maxX, offsetAdjuster)

    if showFeatures:
        if application == 'blastx':
            featureAdder = ProteinFeatureAdder()
        else:
            featureAdder = NucleotideFeatureAdder()

        features = featureAdder.add(featureAx, title, minX, maxX,
                                    offsetAdjuster)

        # If there are features and there weren't too many of them, add
        # vertical feature lines to the reads and ORF axes.
        if features and not featureAdder.tooManyFeaturesToPlot:
            for feature in features:
                start = feature.start()
                end = feature.end()
                color = feature.color
                readsAx.axvline(x=start, color=color)
                readsAx.axvline(x=end, color='#cccccc')
                if showOrfs:
                    orfAx.axvline(x=start, color=color)
                    orfAx.axvline(x=end, color='#cccccc')
                    orfReversedAx.axvline(x=start, color=color)
                    orfReversedAx.axvline(x=end, color='#cccccc')
    else:
        features = None

    # We'll return some information we've gathered.
    result = {
        'features': features,
    }

    if plot == 'e values' and plotInfo['zeroEValueFound']:
        # Add the horizontal divider between the highest e-value and the
        # randomly higher ones (if any).
        readsAx.axhline(y=maxY + 0.5, color='#cccccc', linewidth=0.5)

    # Titles, axis, etc.
    if createFigure:
        readCount = blastHits.titles[title]['readCount']
        hspTotal = plotInfo['hspTotal']
        figure.suptitle(
            '%s\nLength %d %s, %d read%s hit, %d HSP%s in total (%d shown).' %
            (
                sequence.description,
                len(sequence), 'aa' if application == 'blastx' else 'nt',
                readCount, '' if readCount == 1 else 's',
                hspTotal, '' if hspTotal == 1 else 's',
                len(plotInfo['items'])
            ),
            fontsize=20)

    # Add a title and y-axis label, but only if we made the reads axes.
    if createdReadsAx:
        readsAx.set_title('Read alignments', fontsize=20)
        if params['rankValues']:
            if plot == 'bit scores':
                plt.ylabel('bit score rank', fontsize=17)
            else:
                plt.ylabel('e-value rank', fontsize=17)
        else:
            if plot == 'bit scores':
                plt.ylabel('Bit score', fontsize=17)
            else:
                plt.ylabel('$- log_{10}(e)$', fontsize=17)

    # Set the x-axis limits.
    if xRange == 'subject':
        readsAx.set_xlim([minX - 1, maxX + 1])
    else:
        # Look at all the HSPs for this subject and figure out the min read
        # start and max read end so we can set the X-axis.
        first = True
        for item in plotInfo['items']:
            hsp = item['hsp']
            if first:
                queryMin = hsp['queryStart']
                queryMax = hsp['queryEnd']
                first = False
            else:
                if hsp['queryStart'] < queryMin:
                    queryMin = hsp['queryStart']
                if hsp['queryEnd'] > queryMax:
                    queryMax = hsp['queryEnd']

        readsAx.set_xlim([queryMin, queryMax])
        result['queryMin'] = queryMin
        result['queryMax'] = queryMax

    readsAx.set_ylim([0, maxYIncludingRandoms + 1])
    readsAx.grid()
    if createFigure:
        if showFigure:
            plt.show()
        if imageFile:
            figure.savefig(imageFile)
    stop = time()
    if not quiet:
        report('Graph generated in %.3f mins.' % ((stop - start) / 60.0))

    return result


def alignmentPanel(blastHits, sortOn='bitScoreMax', interactive=True,
                   outputDir=None, idList=False, equalizeXAxes=True,
                   xRange='subject', plot='e values'):
    """
    Produces a rectangular panel of graphs that each contain an alignment graph
    against a given sequence.

    @param blastHits: A L{dark.blast.BlastHits} instance.
    sortOn: The attribute to sort subplots on. Either "eMean",
        "eMedian", "eMin", "bitScoreMax", "bitScoreMean", "bitScoreMedian",
        "readCount", "length", or "title".
    interactive: If C{True}, we are interactive and should display the panel
        using figure.show etc.
    outputDir: If not None, specifies a directory to write an HTML summary to.
        If the directory does not exist it will be created.
    idList: a dictionary. The keys is a color and the values is a list of
        read identifiers that should be colored in the respective color.
    equalizeXAxes: if C{True}, adjust the X axis on each alignment plot to be
        the same.
    xRange: set to either 'subject' or 'reads' to indicate the range of the
        X axis.
    @param plot: A C{str}, either 'bit scores' or 'e values', to indicate
        which values should be plotted on the graph Y axis.
    """

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    assert plot in ('bit scores', 'e values'), (
        'plot must be either "bit scores", or "e values".')

    if not (interactive or outputDir):
        raise ValueError('Either interactive or outputDir must be True')

    start = time()
    titles = blastHits.sortTitlesOnPlotInfo(sortOn)
    cols = 5
    rows = int(len(titles) / cols) + (0 if len(titles) % cols == 0 else 1)
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    report('Plotting %d titles in %dx%d grid, sorted on %s' %
           (len(titles), rows, cols, sortOn))

    maxYIncludingRandoms = maxY = minY = maxX = minX = None
    queryMax = queryMin = None

    if outputDir:
        if os.access(outputDir, os.F_OK):
            # outputDir exists. Check it's a directory.
            mode = os.stat(outputDir).st_mode
            assert S_ISDIR(mode), "%r is not a directory." % outputDir
        else:
            os.mkdir(outputDir)
        htmlOutput = AlignmentPanelHTML(outputDir, blastHits)

    coords = dimensionalIterator((rows, cols))

    for i, title in enumerate(titles):
        plotInfo = blastHits.titles[title]['plotInfo']
        row, col = coords.next()
        print '%d: %s %s' % (i, title, NCBISequenceLinkURL(title, ''))
        if interactive:
            alignmentInfo = alignmentGraph(
                blastHits, title, addQueryLines=True, showFeatures=False,
                colorQueryBases=False, createFigure=False, showFigure=False,
                readsAx=ax[row][col], quiet=True, idList=idList, xRange=xRange,
                plot=plot, showOrfs=False)

        if outputDir:
            imageBasename = '%d.png' % i
            imageFile = '%s/%s' % (outputDir, imageBasename)
            alignmentInfo = alignmentGraph(
                blastHits, title, addQueryLines=True, showFeatures=True,
                colorQueryBases=False, showFigure=False, imageFile=imageFile,
                quiet=True, idList=idList, xRange=xRange, plot=plot,
                showOrfs=True)

            # Close the image plot, otherwise it will be displayed when we
            # call plt.show below.
            plt.close()
            htmlOutput.addImage(imageBasename, title, alignmentInfo, plotInfo)

        readCount = blastHits.titles[title]['readCount']
        hspCount = len(plotInfo['items'])
        plotTitle = ('%d: %s\nLength %d, %d read%s hit.\n%d HSP%s (of %d) '
                     'shown.' % (
                         i, title.split(' ', 1)[1][:40],
                         blastHits.titles[title]['length'],
                         readCount, '' if readCount == 1 else 's',
                         hspCount, '' if hspCount == 1 else 's',
                         plotInfo['hspTotal']))

        if plotInfo['items']:
            if plot == 'bit scores':
                if blastHits.plotParams['rankValues']:
                    plotTitle += '\nY axis is ranked bit score'
                else:
                    median = plotInfo['bitScoreMedian']
                    mean = plotInfo['bitScoreMean']
                    ma = plotInfo['bitScoreMax']
                    plotTitle += '\nmax %.2f, median %.2f, mean %.2f' % (
                        ma, median, mean)
            else:
                if blastHits.plotParams['rankValues']:
                    plotTitle += '\nY axis is ranked e-value'
                else:
                    median = plotInfo['originalEMedian']
                    mean = plotInfo['originalEMean']
                    mi = plotInfo['originalEMin']
                    plotTitle += '\nmin %.2f, median %.2f, mean %.2f' % (
                        mi, median, mean)

        ax[row][col].set_title(plotTitle, fontsize=10)

        if i == 0:
            # This is our first title. Unconditionally set our max/min values.
            if plot == 'bit scores':
                maxYIncludingRandoms = maxY = plotInfo['bitScoreMax']
                minY = plotInfo['bitScoreMin']
            else:
                maxYIncludingRandoms = plotInfo['maxEIncludingRandoms']
                maxY = plotInfo['maxE']
                minY = plotInfo['minE']
            maxX = plotInfo['maxX']
            minX = plotInfo['minX']
            if xRange == 'reads':
                queryMin = alignmentInfo['queryMin']
                queryMax = alignmentInfo['queryMax']
        else:
            # Conditionally adjust max/min values.
            if plot == 'bit scores':
                if plotInfo['bitScoreMax'] > maxY:
                    maxY = maxYIncludingRandoms = plotInfo['bitScoreMax']
                if plotInfo['bitScoreMin'] < minY:
                    minY = plotInfo['bitScoreMin']
            else:
                if plotInfo['maxEIncludingRandoms'] > maxYIncludingRandoms:
                    maxYIncludingRandoms = plotInfo['maxEIncludingRandoms']
                if plotInfo['maxE'] > maxY:
                    maxY = plotInfo['maxE']
                if plotInfo['minE'] < minY:
                    minY = plotInfo['minE']
            if plotInfo['maxX'] > maxX:
                maxX = plotInfo['maxX']
            if plotInfo['minX'] < minX:
                minX = plotInfo['minX']
            if xRange == 'reads':
                if alignmentInfo['queryMin'] < queryMin:
                    queryMin = alignmentInfo['queryMin']
                if alignmentInfo['queryMax'] > queryMax:
                    queryMax = alignmentInfo['queryMax']

    # Post-process graphs to adjust axes, etc.

    if xRange == 'reads':
        # We're showing the reads range on the X axis.
        minX, maxX = queryMin, queryMax

    coords = dimensionalIterator((rows, cols))
    for title in titles:
        row, col = coords.next()
        a = ax[row][col]
        a.set_ylim([0, maxYIncludingRandoms + 1])
        if equalizeXAxes:
            a.set_xlim([minX, maxX])
        a.set_yticks([])
        a.set_xticks([])
        plotInfo = blastHits.titles[title]['plotInfo']

        if xRange == 'subject' and minX < 0:
            # Add a vertical line at x=0 so we can see the 'whiskers' of
            # reads that extend to the left of the sequence we're aligning
            # against.
            a.axvline(x=0, color='#cccccc')

        # Add a line on the right of each sub-plot so we can see where the
        # sequence ends (as all panel graphs have the same width and we
        # otherwise couldn't tell).
        sequenceLen = blastHits.titles[title]['length']
        if blastHits.plotParams['logLinearXAxis']:
            offsetAdjuster = plotInfo['offsetAdjuster'].adjustOffset
            sequenceLen = offsetAdjuster(sequenceLen)
        a.axvline(x=sequenceLen, color='#cccccc')

    # Hide the final panel graphs (if any) that have no content. We do this
    # because the panel is a rectangular grid and some of the plots at the
    # end of the last row can be unused.
    for row, col in coords:
        ax[row][col].axis('off')

    # plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.93,
    # wspace=0.1, hspace=None)
    plt.subplots_adjust(hspace=0.4)
    figure.suptitle('X: %d to %d, Y (%s): %d to %d' %
                    (minX, maxX, plot, int(minY), int(maxY)), fontsize=20)
    figure.set_size_inches(5 * cols, 3 * rows, forward=True)
    if outputDir:
        panelFilename = 'alignment-panel.png'
        figure.savefig('%s/%s' % (outputDir, panelFilename))
        htmlOutput.close(panelFilename)
    if interactive:
        figure.show()
    stop = time()
    report('Alignment panel generated in %.3f mins.' % ((stop - start) / 60.0))


def bitScoreAndEValueScatter(blastHits):
    """
    Show a scatter plot of bit scores against e-values for a set of hit
    sequences.

    @param blastHits: A L{dark.blast.BlastHits} instance.
    """
    bitScores = []
    eValues = []

    for title, hitInfo in blastHits.titles.iteritems():
        if hitInfo['plotInfo'] is None:
            continue
        for item in hitInfo['plotInfo']['items']:
            bitScores.append(item['bitScore'])
            eValues.append(item['convertedE'])

    savedFigsize = matplotlib.rcParams['figure.figsize']
    matplotlib.rcParams['figure.figsize'] = [12, 12]
    plt.gray()
    plt.grid(True)
    plt.scatter(bitScores, eValues, marker='.', s=10, color='red')
    plt.xlabel('Bit score', fontsize=17)
    plt.ylabel('$-log_{10}(evalue)$', fontsize=17)
    plt.title('%d bit scores vs E-values' % len(bitScores), fontsize=20)
    plt.show()
    matplotlib.rcParams['figure.figsize'] = savedFigsize


def evalueGraph(records, rows, cols, find=None, titles=True, minHits=1,
                width=5, height=5):
    """
    Produces a rectangular panel of graphs that each show sorted e-values for
    a read. Read hits against a certain strain (see find, below) are
    highlighted.

    find: A function that can be passed a sequence title. If the function
    returns C{True} a (currently) red dot is put into the graph at that point.
    titles: Show read sequence names.
    minHits: only show reads with at least this many hits.
    """
    f, ax = plt.subplots(rows, cols)
    globalMaxE = 0.0
    globalMaxDescriptions = 0
    recordCount = 0
    done = False
    lowHitCount = 0
    for row in xrange(rows):
        if done:
            break
        for col in xrange(cols):
            try:
                record = records.next()
                while len(record.descriptions) < minHits:
                    # print 'rejecting', record.query
                    lowHitCount += 1
                    record = records.next()
            except StopIteration:
                done = True
                break
            else:
                recordCount += 1
                if len(record.descriptions) > globalMaxDescriptions:
                    globalMaxDescriptions = len(record.descriptions)
                evalues = []
                foundx = []
                foundy = []
                for i, desc in enumerate(record.descriptions):
                    # NOTE: We are looping over the descriptions here, not
                    # the multiple HSPs in the alignments. The description
                    # describes only the first (i.e., the best) HSP.
                    e = -1.0 * log10(record.alignments[i].hsps[0].expect)
                    if e < 0:
                        break
                    evalues.append(e)
                    if find and find(desc.title):
                        foundx.append(i)
                        foundy.append(e)
                a = ax[row][col]
                if evalues:
                    maxE = max(evalues)
                    if maxE > globalMaxE:
                        globalMaxE = maxE
                    x = np.arange(0, len(evalues))
                    a.plot(x, evalues)
                if foundx:
                    # a.plot(foundx, foundy, 'ro', markersize=5)
                    a.plot(foundx, foundy, 'ro')
                if titles:
                    a.set_title('%s (%d)' %
                                (record.query, record.query_length),
                                fontsize=10)

    count = 0
    for row in xrange(rows):
        for col in xrange(cols):
            if count == recordCount:
                break
            count += 1
            a = ax[row][col]
            a.axis([0, globalMaxDescriptions, 0, globalMaxE])
            # a.set_yscale('log')
            a.set_yticks([])
            a.set_xticks([])

    plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.93,
                        wspace=0.1, hspace=None)
    f.suptitle('maxHits %d, maxE %f, ignored %d, (minHits %d)' %
               (globalMaxDescriptions, globalMaxE, lowHitCount, minHits))
    f.set_size_inches(width, height, forward=True)
    # f.set_size_inches(10, 10)
    # f.savefig('evalues.png')
    plt.show()


def scatterAlign(seq1, seq2, window=7):
    """
    Visually align two sequences.
    """
    d1 = defaultdict(list)
    d2 = defaultdict(list)
    for (seq, section_dict) in [(seq1, d1), (seq2, d2)]:
        for i in range(len(seq) - window):
            section = seq[i:i + window]
            section_dict[section].append(i)
    matches = set(d1).intersection(d2)
    print '%i unique matches' % len(matches)
    x = []
    y = []
    for section in matches:
        for i in d1[section]:
            for j in d2[section]:
                x.append(i)
                y.append(j)
    # plt.cla()  # clear any prior graph
    plt.gray()
    plt.scatter(x, y)
    plt.xlim(0, len(seq1) - window)
    plt.ylim(0, len(seq2) - window)
    plt.xlabel('length %i bp' % (len(seq1)))
    plt.ylabel('length %i bp' % (len(seq2)))
    plt.title('Dot plot using window size %i\n(allowing no mis-matches)' %
              window)
    plt.show()

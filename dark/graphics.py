import os
from stat import S_ISDIR
from math import ceil
from collections import defaultdict
from time import ctime, time
from math import log10
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib import gridspec
from IPython.display import HTML

from dark import html
from dark.baseimage import BaseImage
from dark.dimension import dimensionalIterator
from dark.html import NCBISequenceLink
from dark.intervals import ReadIntervals
from dark import genbank
from dark import ncbidb
from dark import features

START_CODONS = set(['ATG'])
STOP_CODONS = set(['TAA', 'TAG', 'TGA'])

QUERY_COLORS = {
    'A': (1.0, 0.0, 0.0),  # Red.
    'C': (0.0, 0.0, 1.0),  # Blue.
    'G': (0.0, 1.0, 0.0),  # Green.
    'N': (1.0, 0.0, 1.0),  # Purple.
    'T': (1.0, 0.8, 0.0),  # Orange.
    'gap': (0.2, 0.2, 0.2),  # Almost black.
    'match': (0.9, 0.9, 0.9),  # Almost white.
}

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
            '%3d: count=%4d, len=%7d, median(e)=%20s mean(e)=%20s: %s' %
            (i, hitInfo['readCount'], hitInfo['length'],
             hitInfo['eMedian'], hitInfo['eMean'], link))
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


def alignmentGraph(blastHits, title, addQueryLines=True, showFeatures=True,
                   colorQueryBases=False, createFigure=True, showFigure=True,
                   readsAx=None, imageFile=None, quiet=False, idList=False,
                   xRange='subject'):
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
    """

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    start = time()

    titleId = title.split(' ', 1)[0]
    sequence = ncbidb.getSequence(titleId, blastHits.records.blastDb)

    if createFigure:
        dpi = 80
        # width = max(float(len(sequence)) / float(dpi), 25)
        width = 20
        figure = plt.figure(figsize=(width, 20), dpi=dpi)

    createdReadsAx = readsAx is None

    if showFeatures:
        gbSeq = genbank.getSequence(titleId)
        gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 12])
        featureAx = plt.subplot(gs[0, 0])
        orfAx = plt.subplot(gs[1, 0])
        orfReversedAx = plt.subplot(gs[2, 0])
        readsAx = readsAx or plt.subplot(gs[3, 0])
    else:
        featureEndpoints = []
        readsAx = readsAx or plt.subplot(111)

    params = blastHits.plotParams
    logLinearXAxis = params['logLinearXAxis']
    rankEValues = params['rankEValues']

    plotInfo = blastHits.titles[title]['plotInfo']
    items = plotInfo['items']
    maxEIncludingRandoms = int(ceil(plotInfo['maxEIncludingRandoms']))
    maxE = int(ceil(plotInfo['maxE']))
    minE = int(plotInfo['minE'])
    maxX = plotInfo['maxX']
    minX = plotInfo['minX']

    # Add light grey vertical rectangles to show the logarithmic gaps. Add
    # these first so that reads will be plotted on top of them. Only draw
    # gaps that are more than 20 pixels wide as we could have millions of
    # tiny gaps for a bacteria and drawing them all will be slow and only
    # serves to make the entire background grey.
    if logLinearXAxis and len(plotInfo['offsetAdjuster'].adjustments()) < 100:
        offsetAdjuster = plotInfo['offsetAdjuster'].adjustOffset
        for (intervalType, interval) in plotInfo['readIntervals'].walk():
            if intervalType == ReadIntervals.EMPTY:
                adjustedStart = offsetAdjuster(interval[0])
                adjustedStop = offsetAdjuster(interval[1])
                width = adjustedStop - adjustedStart
                if width >= SMALLEST_LOGGED_GAP_TO_DISPLAY:
                    rect = Rectangle(
                        (adjustedStart, 0), width,
                        maxEIncludingRandoms + 1, color='#f4f4f4')
                    readsAx.add_patch(rect)

    if colorQueryBases:
        # Color each query by its bases.
        xScale = 3
        yScale = 2
        baseImage = BaseImage(
            maxX - minX,
            maxEIncludingRandoms - minE + (1 if rankEValues else 0),
            xScale, yScale)
        for item in items:
            hsp = item['hsp']
            e = item['convertedE'] - minE
            # If the product of the subject and query frame values is +ve,
            # then they're either both +ve or both -ve, so we just use the
            # query as is. Otherwise, we need to reverse complement it.
            if item['frame']['subject'] * item['frame']['query'] > 0:
                query = blastHits.fasta[item['readNum']].seq
            else:
                # One of the subject or query has negative sense.
                query = blastHits.fasta[
                    item['readNum']].reverse_complement().seq
            queryStart = hsp['queryStart']
            # There are 3 parts of the query string we need to display. 1)
            # the left part (if any) before the matched part of the
            # subject. 2) the matched part (which can include gaps in the
            # query and/or subject). 3) the right part (if any) after the
            # matched part.

            # 1. Left part.
            xOffset = queryStart - minX
            queryOffset = 0
            for queryIndex in xrange(hsp['subjectStart'] - queryStart):
                color = QUERY_COLORS[query[queryOffset + queryIndex]]
                baseImage.set(xOffset + queryIndex, e, color)

            # 2. Match part.
            xOffset = hsp['subjectStart'] - minX
            xIndex = 0
            queryOffset = hsp['subjectStart'] - hsp['queryStart']
            origSubject = item['origHsp'].sbjct
            origQuery = item['origHsp'].query
            for matchIndex in xrange(len(origSubject)):
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
                            color = QUERY_COLORS[origQuery[matchIndex]]
                    baseImage.set(xOffset + xIndex, e, color)
                    xIndex += 1

            # 3. Right part.
            xOffset = hsp['subjectEnd'] - minX
            queryOffset = hsp['subjectEnd'] - hsp['queryStart']
            for queryIndex in xrange(hsp['queryEnd'] - hsp['subjectEnd']):
                color = QUERY_COLORS[query[queryOffset + queryIndex]]
                baseImage.set(xOffset + queryIndex, e, color)

                # a gap in the subject was needed to match a different query
                # than the one we're matching against.
                # this throws off the indexing.
                # if this happen, show a gap when coloring in the base of the
                # query as gap.

        readsAx.imshow(baseImage.data, aspect='auto', origin='lower',
                       interpolation='nearest',
                       extent=[minX, maxX, minE, maxEIncludingRandoms])
    else:
        # Add horizontal lines for all the query sequences. These will be the
        # grey 'whiskers' in the plots once we (below) draw the matched part
        # on top of part of them.
        if addQueryLines:
            for item in items:
                e = item['convertedE']
                hsp = item['hsp']
                line = Line2D([hsp['queryStart'], hsp['queryEnd']], [e, e],
                              color='#aaaaaa')
                readsAx.add_line(line)

        # Add the horizontal BLAST alignment lines.
        for item in items:
            e = item['convertedE']
            hsp = item['hsp']
            line = Line2D([hsp['subjectStart'], hsp['subjectEnd']], [e, e],
                          color='blue')
            readsAx.add_line(line)

        # TODO: Barbara, add some kind of comment here please. I think this
        # could be done more efficiently.
        if idList:
            for item in items:
                for key in idList:
                    for ids in idList[key]:
                        if ids == item['query']:
                            e = item['convertedE']
                            hsp = item['hsp']
                            line = Line2D([hsp['subjectStart'],
                                           hsp['subjectEnd']], [e, e],
                                          color=key)
                            readsAx.add_line(line)

    # Add to ORF figures and add vertical lines for the sequence features.
    # The feature and ORF display are linked and they shouldn't be. This
    # code is a horrible mess.
    if showFeatures:
        if logLinearXAxis:
            offsetAdjuster = plotInfo['offsetAdjuster'].adjustOffset
        else:
            offsetAdjuster = lambda x: x

        featureEndpoints = features.addFeatures(featureAx, gbSeq, minX, maxX,
                                                offsetAdjuster)
        if len(featureEndpoints) < 20:
            for fe in featureEndpoints:
                readsAx.axvline(x=fe['start'], color=fe['color'])
                readsAx.axvline(x=fe['end'], color='#cccccc')
            features.addORFs(orfAx, sequence.seq, minX, maxX, featureEndpoints,
                             offsetAdjuster)
        else:
            features.addORFs(orfAx, sequence.seq, minX, maxX, [],
                             offsetAdjuster)

        features.addReversedORFs(orfReversedAx,
                                 sequence.reverse_complement().seq,
                                 minX, maxX, offsetAdjuster)

    # We'll return some information we've gathered.
    result = {
        'features': featureEndpoints,
    }

    # Add the horizontal divider between the highest e-value and the randomly
    # higher ones (if any).
    if plotInfo['zeroEValueFound']:
        readsAx.axhline(y=maxE + 0.5, color='#cccccc', linewidth=1)

    # Titles, axis, etc.
    if createFigure:
        readCount = blastHits.titles[title]['readCount']
        hspTotal = plotInfo['hspTotal']
        figure.suptitle(
            '%s\nLength %d, %d read%s hit, %d HSP%s in total (%d shown).' %
            (
                sequence.description, len(sequence),
                readCount, '' if readCount == 1 else 's',
                hspTotal, '' if hspTotal == 1 else 's',
                len(plotInfo['items'])
            ),
            fontsize=20)

    if createdReadsAx:
        # Only add title and y-axis label if we made the reads axes.
        readsAx.set_title('Read alignments', fontsize=20)
        if rankEValues:
            plt.ylabel('e value rank', fontsize=17)
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

    readsAx.set_ylim([0, maxEIncludingRandoms + 1])
    readsAx.grid()
    if createFigure:
        if showFigure:
            plt.show()
        if imageFile:
            figure.savefig(imageFile)
    stop = time()
    if not quiet:
        report('Graph generated in %.3f mins. Read count: %d. HSP count: %d.' %
               ((stop - start) / 60.0, len(blastHits.fasta), len(items)))

    return result


def alignmentPanel(blastHits, sortOn='eMedian', interactive=True,
                   outputDir=None, idList=False, equalizeXAxes=True,
                   xRange='subject'):
    """
    Produces a rectangular panel of graphs that each contain an alignment graph
    against a given sequence.

    @param blastHits: A L{dark.blast.BlastHits} instance.
    sortOn: The attribute to sort subplots on. Either "eMean", "eMedian",
        "title" or "reads"
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
    """

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    if not (interactive or outputDir):
        raise ValueError('Either interactive or outputDir must be True')

    start = time()
    titles = blastHits.sortTitles(sortOn)
    cols = 5
    rows = int(len(titles) / cols) + (0 if len(titles) % cols == 0 else 1)
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    report('Plotting %d titles in %dx%d grid, sorted on %s' %
           (len(titles), rows, cols, sortOn))

    maxEIncludingRandoms = maxE = minE = maxX = minX = None
    queryMax = queryMin = None

    if outputDir:
        if os.access(outputDir, os.F_OK):
            # outputDir exists. Check it's a directory.
            mode = os.stat(outputDir).st_mode
            assert S_ISDIR(mode), "%r is not a directory." % outputDir
        else:
            os.mkdir(outputDir)
        htmlOutput = html.AlignmentPanelHTML(outputDir, blastHits)

    coords = dimensionalIterator((rows, cols))

    for i, title in enumerate(titles):
        plotInfo = blastHits.titles[title]['plotInfo']
        row, col = coords.next()
        print '%d: %s %s' % (i, title, html.NCBISequenceLinkURL(title, ''))
        if interactive:
            alignmentInfo = alignmentGraph(
                blastHits, title, addQueryLines=True, showFeatures=False,
                colorQueryBases=False, createFigure=False, showFigure=False,
                readsAx=ax[row][col], quiet=True, idList=idList, xRange=xRange)

        if outputDir:
            imageBasename = '%d.png' % i
            imageFile = '%s/%s' % (outputDir, imageBasename)
            alignmentInfo = alignmentGraph(
                blastHits, title, addQueryLines=True, showFeatures=False,
                colorQueryBases=False, showFigure=False, imageFile=imageFile,
                quiet=True, idList=idList, xRange=xRange)

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
            if blastHits.plotParams['rankEValues']:
                plotTitle += '\ne-values are ranks'
            else:
                eValues = [item['convertedE'] for item in plotInfo['items']]
                median = np.median(eValues)
                mean = np.mean(eValues)
                plotTitle += '\nmedian 1e-%d, mean 1e-%d' % (median, mean)

        ax[row][col].set_title(plotTitle, fontsize=10)

        if i == 0:
            # This is our first title. Unconditionally set our max/min values.
            maxEIncludingRandoms = plotInfo['maxEIncludingRandoms']
            maxE = plotInfo['maxE']
            minE = plotInfo['minE']
            maxX = plotInfo['maxX']
            minX = plotInfo['minX']
            if xRange == 'reads':
                queryMin = alignmentInfo['queryMin']
                queryMax = alignmentInfo['queryMax']
        else:
            # Conditionally adjust max/min values.
            if plotInfo['maxEIncludingRandoms'] > maxEIncludingRandoms:
                maxEIncludingRandoms = plotInfo['maxEIncludingRandoms']
            if plotInfo['maxE'] > maxE:
                maxE = plotInfo['maxE']
            if plotInfo['minE'] < minE:
                minE = plotInfo['minE']
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
        a.set_ylim([0, maxEIncludingRandoms + 1])
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

        # Add light grey vertical bars to show the logarithmic gaps. Add
        # these only in the region above the highest read in the individual
        # plot. If we simply added the bar to the full height of the plot
        # it would obscure the reads it overlapped.
        if (blastHits.plotParams['logLinearXAxis'] and
                len(plotInfo['offsetAdjuster'].adjustments()) < 100):
            thisMaxEIncludingRandoms = plotInfo['maxEIncludingRandoms']
            for (intervalType, interval) in plotInfo['readIntervals'].walk():
                if intervalType == ReadIntervals.EMPTY:
                    adjustedStart = offsetAdjuster(interval[0])
                    adjustedStop = offsetAdjuster(interval[1])
                    width = adjustedStop - adjustedStart
                    if width >= SMALLEST_LOGGED_GAP_TO_DISPLAY:
                        height = (maxEIncludingRandoms -
                                  thisMaxEIncludingRandoms)
                        rect = Rectangle(
                            (adjustedStart, thisMaxEIncludingRandoms + 1),
                            width, height, color='#f4f4f4')
                        a.add_patch(rect)

    # Hide the final panel graphs (if any) that have no content. We do this
    # because the panel is a rectangular grid and some of the plots at the
    # end of the last row can be unused.
    for row, col in coords:
        ax[row][col].axis('off')

    # plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.93,
    # wspace=0.1, hspace=None)
    plt.subplots_adjust(hspace=0.4)
    figure.suptitle('X: %d to %d, Y: %d to %d' %
                    (minX, maxX, int(minE), int(maxE)), fontsize=20)
    figure.set_size_inches(5 * cols, 3 * rows, forward=True)
    if outputDir:
        panelFilename = 'alignment-panel.png'
        figure.savefig('%s/%s' % (outputDir, panelFilename))
        htmlOutput.close(panelFilename)
    if interactive:
        figure.show()
    stop = time()
    report('Alignment panel generated in %.3f mins.' % ((stop - start) / 60.0))


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
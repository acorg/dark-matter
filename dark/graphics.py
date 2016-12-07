import os
from copy import deepcopy
from stat import S_ISDIR
from math import ceil
from collections import defaultdict
from time import ctime, time
from textwrap import fill

try:
    import matplotlib.pyplot as plt
except ImportError:
    import platform
    if platform.python_implementation() == 'PyPy':
        # PyPy doesn't have a version of matplotlib. Make fake classes and
        # a Line2D function that raises if used. This allows us to use
        # other 'dark' code that happens to import dark.graphics but which
        # does not use the functions that rely on matplotlib.
        class plt(object):
            def __getattr__(self, _):
                raise NotImplementedError(
                    'matplotlib is not supported under pypy')

        gridspec = patches = plt

        def Line2D(*args, **kwargs):
            raise NotImplementedError('matplotlib is not supported under pypy')
    else:
        raise
else:
    from matplotlib.lines import Line2D
    from matplotlib import gridspec, patches

import numpy as np

from dark.aa import propertiesForSequence, clustersForSequence
from dark.baseimage import BaseImage
from dark.dimension import dimensionalIterator
from dark.html import AlignmentPanelHTMLWriter, NCBISequenceLinkURL
from dark.intervals import ReadIntervals
from dark.features import ProteinFeatureAdder, NucleotideFeatureAdder
from dark import orfs
from dark.intervals import OffsetAdjuster
from dark.score import HigherIsBetterScore


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

# From http://www.randalolson.com/2014/06/28/how-to-make-beautiful-data-\
# visualizations-in-python-with-matplotlib/
#
# These are the "Tableau 20" colors as RGB.
TABLEAU20 = [
    (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
    (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
    (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
    (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
    (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the above RGB values to the [0, 1] range, the format matplotlib
# accepts.
for i in range(len(TABLEAU20)):
    r, g, b = TABLEAU20[i]
    TABLEAU20[i] = (r / 255.0, g / 255.0, b / 255.0)

# If we're making a plot that has a log-linear X axis, don't show
# background light grey rectangles for any gap whose (logged) width is less
# than SMALLEST_LOGGED_GAP_TO_DISPLAY.
SMALLEST_LOGGED_GAP_TO_DISPLAY = 20

# Y (Score) axis extra spacing above the best read match. The score of the
# best read HSP will be multiplied by this value, the result converted to
# an int, and then used as the upper bound on the Y axis. Adding 1% seems
# to be a good heuristic.
Y_AXIS_UPPER_PADDING = 1.01

# The default base of the logarithm to use when logLinearXAxis is used to
# produce an alignment graph.
DEFAULT_LOG_LINEAR_X_AXIS_BASE = 1.1


def report(msg):
    print('%s: %s' % (ctime(time()), msg))


def alignmentGraph(titlesAlignments, title, addQueryLines=True,
                   showFeatures=True, logLinearXAxis=False,
                   logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE, rankScores=False,
                   colorQueryBases=False, createFigure=True, showFigure=True,
                   readsAx=None, imageFile=None, quiet=False, idList=False,
                   xRange='subject', showOrfs=True):
    """
    Align a set of matching reads against a BLAST hit.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param title: A C{str} sequence title that was hit by BLAST. We plot the
        reads that hit this title.
    @param addQueryLines: if C{True}, draw query lines in full (these will then
        be partly overdrawn by the HSP match against the subject). These are
        the 'whiskers' that potentially protrude from each side of a query.
    @param showFeatures: if C{True}, look online for features of the subject
        sequence (given by hitId).
    @param logLinearXAxis: if C{True}, convert read offsets so that empty
        regions in the plot we're preparing will only be as wide as their
        logged actual values.
    @param logBase: The base of the logarithm to use if logLinearXAxis is
        C{True}.
    @param: rankScores: If C{True}, change the e-values and bit scores for the
        reads for each title to be their rank (worst to best).
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
    @param showOrfs: If C{True}, open reading frames will be displayed.
    """

    startTime = time()

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    if createFigure:
        width = 20
        figure = plt.figure(figsize=(width, 20))

    createdReadsAx = readsAx is None

    if showFeatures:
        if showOrfs:
            gs = gridspec.GridSpec(4, 1, height_ratios=[3, 1, 1, 12])
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

    # Make a deep copy of the title alignments. We're potentially going to
    # change the HSP scores, the X axis offsets, etc., and we don't want to
    # interfere with the data we were passed.
    titleAlignments = deepcopy(titlesAlignments[title])

    readsAlignments = titlesAlignments.readsAlignments
    subjectIsNucleotides = readsAlignments.params.subjectIsNucleotides
    subject = readsAlignments.getSubjectSequence(title)

    if showOrfs and not subjectIsNucleotides:
        # We cannot show ORFs when displaying protein plots.
        showOrfs = False

    # Allow the class of titlesAlignments to adjust HSPs for plotting,
    # if it has a method for doing so.
    try:
        adjuster = readsAlignments.adjustHspsForPlotting
    except AttributeError:
        pass
    else:
        adjuster(titleAlignments)

    if rankScores:
        reverse = titlesAlignments.scoreClass is not HigherIsBetterScore
        for rank, hsp in enumerate(sorted(titleAlignments.hsps(),
                                   reverse=reverse), start=1):
            hsp.score.score = rank

    if logLinearXAxis:
        readIntervals = ReadIntervals(titleAlignments.subjectLength)
        # Examine all HSPs so we can build an offset adjuster.
        for hsp in titleAlignments.hsps():
            readIntervals.add(hsp.readStartInSubject, hsp.readEndInSubject)
        # Now adjust offsets in all HSPs.
        offsetAdjuster = OffsetAdjuster(readIntervals, base=logBase)
        for hsp in titleAlignments.hsps():
            offsetAdjuster.adjustHSP(hsp)
        # A function for adjusting other offsets, below.
        adjustOffset = offsetAdjuster.adjustOffset
    else:
        adjustOffset = lambda x: x

    # It would be more efficient to only walk through all HSPs once and
    # compute these values all at once, but for now this is simple and clear.
    maxY = int(ceil(titleAlignments.bestHsp().score.score))
    minY = int(titleAlignments.worstHsp().score.score)
    maxX = max(hsp.readEndInSubject for hsp in titleAlignments.hsps())
    minX = min(hsp.readStartInSubject for hsp in titleAlignments.hsps())

    if xRange == 'subject':
        # We'll display a graph for the full subject range. Adjust X axis
        # min/max to make sure we cover at least zero to the sequence length.
        maxX = max(len(subject), maxX)
        minX = min(0, minX)

    # Swap min & max Y values, if needed, as it's possible we are dealing
    # with LSPs but that the score adjuster made numerically greater values
    # for those that were small.
    if maxY < minY:
        (maxY, minY) = (minY, maxY)

    if logLinearXAxis:
        # Adjust minX and maxX if we have gaps at the subject start or end.
        gaps = list(readIntervals.walk())
        if gaps:
            # Check start of first gap:
            intervalType, (start, stop) = gaps[0]
            if intervalType == ReadIntervals.EMPTY:
                adjustedStart = adjustOffset(start)
                if adjustedStart < minX:
                    minX = adjustedStart
            # Check stop of last gap:
            intervalType, (start, stop) = gaps[-1]
            if intervalType == ReadIntervals.EMPTY:
                adjustedStop = adjustOffset(stop)
                if adjustedStop > maxX:
                    maxX = adjustedStop

    # We're all set up to start plotting the graph.

    # Add light grey vertical rectangles to show the logarithmic gaps. Add
    # these first so that reads will be plotted on top of them. Only draw
    # gaps that are more than SMALLEST_LOGGED_GAP_TO_DISPLAY pixels wide as
    # we could have millions of tiny gaps for a bacteria and drawing them
    # all will be slow and only serves to make the entire background grey.
    if logLinearXAxis and len(offsetAdjuster.adjustments()) < 100:
        for (intervalType, interval) in readIntervals.walk():
            if intervalType == ReadIntervals.EMPTY:
                adjustedStart = adjustOffset(interval[0])
                adjustedStop = adjustOffset(interval[1])
                width = adjustedStop - adjustedStart
                if width >= SMALLEST_LOGGED_GAP_TO_DISPLAY:
                    readsAx.axvspan(adjustedStart, adjustedStop,
                                    color='#f4f4f4')

    if colorQueryBases:
        # Color each query by its bases.
        xScale = 3
        yScale = 2
        baseImage = BaseImage(
            maxX - minX, maxY - minY + (1 if rankScores else 0),
            xScale, yScale)
        for alignment in titleAlignments:
            for hsp in alignment.hsps:
                y = hsp.score.score - minY
                # If the product of the subject and read frame values is +ve,
                # then they're either both +ve or both -ve, so we just use the
                # read as is. Otherwise, we need to reverse complement it.
                if hsp.subjectFrame * hsp.readFrame > 0:
                    query = alignment.read.sequence
                else:
                    # One of the subject or query has negative sense.
                    query = alignment.read.reverseComplement().sequence
                readStartInSubject = hsp.readStartInSubject
                # There are 3 parts of the query string we need to
                # display. 1) the left part (if any) before the matched
                # part of the subject.  2) the matched part (which can
                # include gaps in the query and/or subject). 3) the right
                # part (if any) after the matched part.  For each part,
                # calculate the ranges in which we have to make the
                # comparison between subject and query.

                # NOTE: never use hsp['origHsp'].gaps to calculate the number
                # of gaps, as this number contains gaps in both subject and
                # query.

                # 1. Left part:
                leftRange = hsp.subjectStart - readStartInSubject

                # 2. Match, middle part:
                middleRange = len(hsp.readMatchedSequence)

                # 3. Right part:
                # Using hsp.readEndInSubject - hsp.subjectEnd to calculate the
                # length of the right part leads to the part being too long.
                # The number of gaps needs to be subtracted to get the right
                # length.
                origQuery = hsp.readMatchedSequence.upper()
                rightRange = (hsp.readEndInSubject - hsp.subjectEnd -
                              origQuery.count('-'))

                # 1. Left part.
                xOffset = readStartInSubject - minX
                queryOffset = 0
                for queryIndex in range(leftRange):
                    color = QUERY_COLORS.get(query[queryOffset + queryIndex],
                                             DEFAULT_BASE_COLOR)
                    baseImage.set(xOffset + queryIndex, y, color)

                # 2. Match part.
                xOffset = hsp.subjectStart - minX
                xIndex = 0
                queryOffset = hsp.subjectStart - hsp.readStartInSubject
                origSubject = hsp.subjectMatchedSequence
                for matchIndex in range(middleRange):
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
                                # Query doesn't match subject (and is not a
                                # gap).
                                color = QUERY_COLORS.get(origQuery[matchIndex],
                                                         DEFAULT_BASE_COLOR)
                        baseImage.set(xOffset + xIndex, y, color)
                        xIndex += 1

                # 3. Right part.
                xOffset = hsp.subjectEnd - minX
                backQuery = query[-rightRange:].upper()
                for queryIndex in range(rightRange):
                    color = QUERY_COLORS.get(backQuery[queryIndex],
                                             DEFAULT_BASE_COLOR)
                    baseImage.set(xOffset + queryIndex, y, color)

        readsAx.imshow(baseImage.data, aspect='auto', origin='lower',
                       interpolation='nearest',
                       extent=[minX, maxX, minY, maxY])
    else:
        # Add horizontal lines for all the query sequences. These will be the
        # grey 'whiskers' in the plots once we (below) draw the matched part
        # on top of part of them.
        if addQueryLines:
            for hsp in titleAlignments.hsps():
                y = hsp.score.score
                line = Line2D([hsp.readStartInSubject, hsp.readEndInSubject],
                              [y, y], color='#aaaaaa')
                readsAx.add_line(line)

        # Add the horizontal BLAST alignment lines.

        # If an idList is given set things up to look up read colors.
        readColor = {}
        if idList:
            for color, reads in idList.items():
                for read in reads:
                    if read in readColor:
                        raise ValueError('Read %s is specified multiple '
                                         'times in idList' % read)
                    else:
                        readColor[read] = color

        # Draw the matched region.
        for titleAlignment in titleAlignments:
            readId = titleAlignment.read.id
            for hsp in titleAlignment.hsps:
                y = hsp.score.score
                line = Line2D([hsp.subjectStart, hsp.subjectEnd], [y, y],
                              color=readColor.get(readId, 'blue'))
                readsAx.add_line(line)

    if showOrfs:
        orfs.addORFs(orfAx, subject.sequence, minX, maxX, adjustOffset)
        orfs.addReversedORFs(orfReversedAx,
                             subject.reverseComplement().sequence,
                             minX, maxX, adjustOffset)

    if showFeatures:
        if subjectIsNucleotides:
            featureAdder = NucleotideFeatureAdder()
        else:
            featureAdder = ProteinFeatureAdder()

        features = featureAdder.add(featureAx, title, minX, maxX,
                                    adjustOffset)

        # If there are features and there weren't too many of them, add
        # vertical feature lines to the reads and ORF axes.
        if features and not featureAdder.tooManyFeaturesToPlot:
            for feature in features:
                start = feature.start
                end = feature.end
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
        'adjustOffset': adjustOffset,
        'features': features,
        'minX': minX,
        'maxX': maxX,
        'minY': minY,
        'maxY': maxY,
    }

    # Allow the class of titlesAlignments to add to the plot, if it has a
    # method for doing so.
    try:
        adjuster = readsAlignments.adjustPlot
    except AttributeError:
        pass
    else:
        adjuster(readsAx)

    # Titles, axis, etc.
    if createFigure:
        readCount = titleAlignments.readCount()
        hspCount = titleAlignments.hspCount()
        figure.suptitle(
            '%s\nLength %d %s, %d read%s, %d HSP%s.' %
            (
                fill(subject.id, 80),
                len(subject), 'nt' if subjectIsNucleotides else 'aa',
                readCount, '' if readCount == 1 else 's',
                hspCount, '' if hspCount == 1 else 's'
            ),
            fontsize=20)

    # Add a title and y-axis label, but only if we made the reads axes.
    if createdReadsAx:
        readsAx.set_title('Read alignments', fontsize=20)
        ylabel = readsAlignments.params.scoreTitle
        if rankScores:
            ylabel += ' rank'
        plt.ylabel(ylabel, fontsize=17)

    # Set the x-axis limits.
    readsAx.set_xlim([minX - 1, maxX + 1])

    readsAx.set_ylim([0, int(maxY * Y_AXIS_UPPER_PADDING)])
    readsAx.grid()
    if createFigure:
        if showFigure:
            plt.show()
        if imageFile:
            figure.savefig(imageFile)
    stop = time()
    if not quiet:
        report('Graph generated in %.3f mins.' % ((stop - startTime) / 60.0))

    return result


def alignmentPanel(titlesAlignments, sortOn='maxScore', idList=False,
                   equalizeXAxes=False, xRange='subject', logLinearXAxis=False,
                   rankScores=False, showFeatures=True,
                   logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE):
    """
    Produces a rectangular panel of graphs that each contain an alignment graph
    against a given sequence.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param sortOn: The attribute to sort subplots on. Either "maxScore",
        "medianScore", "readCount", "length", or "title".
    @param idList: A dictionary. Keys are colors and values are lists of read
        ids that should be colored using that color.
    @param equalizeXAxes: If C{True}, adjust the X axis on each alignment plot
        to be the same.
    @param xRange: Set to either 'subject' or 'reads' to indicate the range of
        the X axis.
    @param logLinearXAxis: If C{True}, convert read offsets so that empty
        regions in the plots we're preparing will only be as wide as their
        logged actual values.
    @param logBase: The logarithm base to use if logLinearXAxis is C{True}.
    @param: rankScores: If C{True}, change the scores for the reads for each
        title to be their rank (worst to best).
    @param showFeatures: If C{True}, look online for features of the subject
        sequences.
    @raise ValueError: If C{outputDir} exists but is not a directory or if
        C{xRange} is not "subject" or "reads".
    """

    if xRange not in ('subject', 'reads'):
        raise ValueError('xRange must be either "subject" or "reads".')

    start = time()
    titles = titlesAlignments.sortTitles(sortOn)
    cols = 5
    rows = int(len(titles) / cols) + (0 if len(titles) % cols == 0 else 1)
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    allGraphInfo = {}
    coords = dimensionalIterator((rows, cols))

    report('Plotting %d titles in %dx%d grid, sorted on %s' %
           (len(titles), rows, cols, sortOn))

    for i, title in enumerate(titles):
        titleAlignments = titlesAlignments[title]
        row, col = next(coords)
        report('%d: %s %s' % (i, title, NCBISequenceLinkURL(title, '')))

        # Add a small plot to the alignment panel.
        graphInfo = alignmentGraph(
            titlesAlignments, title, addQueryLines=True,
            showFeatures=showFeatures, rankScores=rankScores,
            logLinearXAxis=logLinearXAxis, logBase=logBase,
            colorQueryBases=False, createFigure=False, showFigure=False,
            readsAx=ax[row][col], quiet=True, idList=idList, xRange=xRange,
            showOrfs=False)

        allGraphInfo[title] = graphInfo
        readCount = titleAlignments.readCount()
        hspCount = titleAlignments.hspCount()

        # Make a short title for the small panel blue plot, ignoring any
        # leading NCBI gi / accession numbers.
        if title.startswith('gi|') and title.find(' ') > -1:
            shortTitle = title.split(' ', 1)[1][:40]
        else:
            shortTitle = title[:40]

        plotTitle = ('%d: %s\nLength %d, %d read%s, %d HSP%s.' % (
            i, shortTitle, titleAlignments.subjectLength,
            readCount, '' if readCount == 1 else 's',
            hspCount, '' if hspCount == 1 else 's'))

        if hspCount:
            if rankScores:
                plotTitle += '\nY axis is ranked score'
            else:
                plotTitle += '\nmax %.2f, median %.2f' % (
                    titleAlignments.bestHsp().score.score,
                    titleAlignments.medianScore())

        ax[row][col].set_title(plotTitle, fontsize=10)

    maxX = max(graphInfo['maxX'] for graphInfo in allGraphInfo.values())
    minX = min(graphInfo['minX'] for graphInfo in allGraphInfo.values())
    maxY = max(graphInfo['maxY'] for graphInfo in allGraphInfo.values())
    minY = min(graphInfo['minY'] for graphInfo in allGraphInfo.values())

    # Post-process graphs to adjust axes, etc.

    coords = dimensionalIterator((rows, cols))
    for title in titles:
        titleAlignments = titlesAlignments[title]
        row, col = next(coords)
        a = ax[row][col]
        a.set_ylim([0, int(maxY * Y_AXIS_UPPER_PADDING)])
        if equalizeXAxes:
            a.set_xlim([minX, maxX])
        a.set_yticks([])
        a.set_xticks([])

        if xRange == 'subject' and minX < 0:
            # Add a vertical line at x=0 so we can see the 'whiskers' of
            # reads that extend to the left of the sequence we're aligning
            # against.
            a.axvline(x=0, color='#cccccc')

        # Add a line on the right of each sub-plot so we can see where the
        # sequence ends (as all panel graphs have the same width and we
        # otherwise couldn't tell).
        sequenceLen = titleAlignments.subjectLength
        if logLinearXAxis:
            sequenceLen = allGraphInfo[title]['adjustOffset'](sequenceLen)
        a.axvline(x=sequenceLen, color='#cccccc')

    # Hide the final panel graphs (if any) that have no content. We do this
    # because the panel is a rectangular grid and some of the plots at the
    # end of the last row may be unused.
    for row, col in coords:
        ax[row][col].axis('off')

    # plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.93,
    # wspace=0.1, hspace=None)
    plt.subplots_adjust(hspace=0.4)
    figure.suptitle('X: %d to %d, Y (%s): %d to %d' %
                    (minX, maxX,
                     titlesAlignments.readsAlignments.params.scoreTitle,
                     int(minY), int(maxY)), fontsize=20)
    figure.set_size_inches(5 * cols, 3 * rows, forward=True)
    figure.show()
    stop = time()
    report('Alignment panel generated in %.3f mins.' % ((stop - start) / 60.0))


def alignmentPanelHTML(titlesAlignments, sortOn='maxScore',
                       outputDir=None, idList=False, equalizeXAxes=False,
                       xRange='subject', logLinearXAxis=False,
                       logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                       rankScores=False, showFeatures=True):
    """
    Produces an HTML index file in C{outputDir} and a collection of alignment
    graphs and FASTA files to summarize the information in C{titlesAlignments}.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param sortOn: The attribute to sort subplots on. Either "maxScore",
        "medianScore", "readCount", "length", or "title".
    @param outputDir: Specifies a C{str} directory to write the HTML to. If
        the directory does not exist it will be created.
    @param idList: A dictionary. Keys are colors and values are lists of read
        ids that should be colored using that color.
    @param equalizeXAxes: If C{True}, adjust the X axis on each alignment plot
        to be the same.
    @param xRange: Set to either 'subject' or 'reads' to indicate the range of
        the X axis.
    @param logLinearXAxis: If C{True}, convert read offsets so that empty
        regions in the plots we're preparing will only be as wide as their
        logged actual values.
    @param logBase: The logarithm base to use if logLinearXAxis is C{True}.
    @param: rankScores: If C{True}, change the scores for the reads for each
        title to be their rank (worst to best).
    @param showFeatures: If C{True}, look online for features of the subject
        sequences.
    @raise TypeError: If C{outputDir} is C{None}.
    @raise ValueError: If C{outputDir} exists but is not a directory or if
        C{xRange} is not "subject" or "reads".
    """

    if xRange not in ('subject', 'reads'):
        raise ValueError('xRange must be either "subject" or "reads".')

    if equalizeXAxes:
        raise NotImplementedError('This feature is not yet implemented.')

    titles = titlesAlignments.sortTitles(sortOn)

    if os.access(outputDir, os.F_OK):
        # outputDir exists. Check it's a directory.
        if not S_ISDIR(os.stat(outputDir).st_mode):
            raise ValueError("%r is not a directory." % outputDir)
    else:
        os.mkdir(outputDir)

    htmlWriter = AlignmentPanelHTMLWriter(outputDir, titlesAlignments)

    for i, title in enumerate(titles):
        # titleAlignments = titlesAlignments[title]

        # If we are writing data to a file too, create a separate file with
        # a plot (this will be linked from the summary HTML).
        imageBasename = '%d.png' % i
        imageFile = '%s/%s' % (outputDir, imageBasename)
        graphInfo = alignmentGraph(
            titlesAlignments, title, addQueryLines=True,
            showFeatures=showFeatures, rankScores=rankScores,
            logLinearXAxis=logLinearXAxis, logBase=logBase,
            colorQueryBases=False, showFigure=False, imageFile=imageFile,
            quiet=True, idList=idList, xRange=xRange, showOrfs=True)

        # Close the image plot to make sure memory is flushed.
        plt.close()
        htmlWriter.addImage(imageBasename, title, graphInfo)

    htmlWriter.close()


def scoreGraph(titlesAlignments, find=None, showTitles=False, figureWidth=5,
               figureHeight=5):
    """
    NOTE: This function has probably bit rotted (but only a little).

    Produce a rectangular panel of graphs, each of which shows sorted scores
    for a title. Matches against a certain sequence title, as determined by
    C{find}, (see below) are highlighted.

    @param find: A function that can be passed a sequence title. If the
        function returns C{True} a red dot is put into the graph at that point
        to highlight the match.
    @param showTitles: If C{True} display read sequence names. The panel tends
        to look terrible when titles are displayed. If C{False}, show no title.
    @param figureWidth: The C{float} width of the figure, in inches.
    @param figureHeight: The C{float} height of the figure, in inches.
    """
    maxScore = None
    maxHsps = 0
    cols = 5
    rows = int(len(titlesAlignments) / cols) + (
        0 if len(titlesAlignments) % cols == 0 else 1)
    f, ax = plt.subplots(rows, cols)
    coords = dimensionalIterator((rows, cols))

    for title in titlesAlignments:
        titleAlignments = titlesAlignments[title]
        row, col = next(coords)
        hspCount = titleAlignments.hspCount()
        if hspCount > maxHsps:
            maxHsps = hspCount
        scores = []
        highlightX = []
        highlightY = []
        for x, titleAlignment in enumerate(titleAlignments):
            score = titleAlignment.hsps[0].score.score
            scores.append(score)
            if find and find(titleAlignment.subjectTitle):
                highlightX.append(x)
                highlightY.append(score)
        a = ax[row][col]
        if scores:
            max_ = max(scores)
            if maxScore is None or max_ > maxScore:
                maxScore = max_
            x = np.arange(0, len(scores))
            a.plot(x, scores)
        if highlightX:
            a.plot(highlightX, highlightY, 'ro')
        if showTitles:
            a.set_title('%s' % title, fontsize=10)

    # Adjust all plots to have the same dimensions.
    coords = dimensionalIterator((rows, cols))
    for _ in range(len(titlesAlignments)):
        row, col = next(coords)
        a = ax[row][col]
        a.axis([0, maxHsps, 0, maxScore])
        # a.set_yscale('log')
        a.set_yticks([])
        a.set_xticks([])

    # Hide the final panel graphs (if any) that have no content. We do this
    # because the panel is a rectangular grid and some of the plots at the
    # end of the last row may be unused.
    for row, col in coords:
        ax[row][col].axis('off')

    plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.93,
                        wspace=0.1, hspace=None)
    f.suptitle('max HSPs %d, max score %f' % (maxHsps, maxScore))
    f.set_size_inches(figureWidth, figureHeight, forward=True)
    # f.savefig('scores.png')
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
    print('%i unique matches' % len(matches))
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


def plotAAProperties(sequence, propertyNames, showLines=True, showFigure=True):
    """
    Plot amino acid property values for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{dark.aa.PROPERTY_DETAILS} C{dict}).
    @param showLines: If C{True}, lines will be drawn between successive AA
        property values. If not, just the values will be plotted as a scatter
        plot (this greatly reduces visual clutter if the sequence is long and
        AA property values are variable).
    @param showFigure: If C{True}, display the plot. Passing C{False} is useful
        in testing.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: The return value from calling dark.aa.propertiesForSequence:
        a C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding property value according to sequence
        position.
    """
    MISSING_AA_VALUE = -1.1
    propertyValues = propertiesForSequence(sequence, propertyNames,
                                           missingAAValue=MISSING_AA_VALUE)
    if showFigure:
        legend = []
        x = np.arange(0, len(sequence))
        plot = plt.plot if showLines else plt.scatter

        for index, propertyName in enumerate(propertyValues):
            color = TABLEAU20[index]
            plot(x, propertyValues[propertyName], color=color)
            legend.append(patches.Patch(color=color, label=propertyName))

        plt.legend(handles=legend, loc=(0, 1.1))
        plt.xlim(-0.2, len(sequence) - 0.8)
        plt.ylim(min(MISSING_AA_VALUE, -1.1), 1.1)
        plt.xlabel('Sequence index')
        plt.ylabel('Property value')
        plt.title(sequence.id)
        plt.show()

    return propertyValues


def plotAAClusters(sequence, propertyNames, showLines=True, showFigure=True):
    """
    Plot amino acid property cluster numbers for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{dark.aa.PROPERTY_CLUSTERS} C{dict}).
    @param showLines: If C{True}, lines will be drawn between successive AA
        property values. If not, just the values will be plotted as a scatter
        plot (this greatly reduces visual clutter if the sequence is long and
        AA property values are variable).
    @param showFigure: If C{True}, display the plot. Passing C{False} is useful
        in testing.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: The return value from calling dark.aa.clustersForSequence:
        a C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding property value according to sequence
        position.
    """
    MISSING_AA_VALUE = 0
    propertyClusters = clustersForSequence(sequence, propertyNames,
                                           missingAAValue=MISSING_AA_VALUE)
    if showFigure:
        minCluster = 1
        maxCluster = -1
        legend = []
        x = np.arange(0, len(sequence))
        plot = plt.plot if showLines else plt.scatter

        for index, propertyName in enumerate(propertyClusters):
            color = TABLEAU20[index]
            clusterNumbers = propertyClusters[propertyName]
            plot(x, clusterNumbers, color=color)
            legend.append(patches.Patch(color=color, label=propertyName))

            propertyMinCluster = min(clusterNumbers)
            if propertyMinCluster < minCluster:
                minCluster = propertyMinCluster

            propertyMaxCluster = max(clusterNumbers)
            if propertyMaxCluster > maxCluster:
                maxCluster = propertyMaxCluster

        plt.legend(handles=legend, loc=(0, 1.1))
        plt.xlim(-0.2, len(sequence) - 0.8)
        plt.ylim(minCluster - 0.5, maxCluster + 0.5)
        plt.yticks(range(maxCluster + 1))
        plt.xlabel('Sequence index')
        plt.ylabel('Property cluster number')
        plt.title(sequence.id)
        plt.show()

    return propertyClusters

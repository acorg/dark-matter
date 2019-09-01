import os
from copy import deepcopy
from stat import S_ISDIR
from math import ceil
from time import ctime, time
from textwrap import fill
from os.path import join

try:
    import matplotlib
    if not os.environ.get('DISPLAY'):
        # Use non-interactive Agg backend
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    import platform
    if platform.python_implementation() == 'PyPy':
        # PyPy doesn't have a version of matplotlib. Make fake classes and
        # a Line2D function that raises if used. This allows us to use
        # other 'dark' code that happens to import dark.civ.graphics but
        # which does not use the functions that rely on matplotlib.
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
    from matplotlib import gridspec

from dark.dimension import dimensionalIterator
from dark.html import NCBISequenceLinkURL
from dark.civ.html import AlignmentPanelHTMLWriter
from dark.intervals import ReadIntervals
from dark.features import ProteinFeatureAdder, NucleotideFeatureAdder
from dark.intervals import OffsetAdjuster
from dark.score import HigherIsBetterScore


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


def alignmentGraph(titlesAlignments, title, accession, addQueryLines=True,
                   showFeatures=True, logLinearXAxis=False,
                   logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE, rankScores=False,
                   createFigure=True, showFigure=True,
                   readsAx=None, imageFile=None, quiet=False, idList=False,
                   xRange='subject'):
    """
    Align a set of matching reads against a BLAST or DIAMOND hit.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param title: A C{str} sequence title that was matched. We plot the
        reads that hit this title.
    @param accession: The C{str} accession number of the matched title.
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

    startTime = time()

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    if createFigure:
        width = 20
        figure = plt.figure(figsize=(width, 20))

    createdReadsAx = readsAx is None

    if showFeatures:
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
        featureAx = plt.subplot(gs[0, 0])
        readsAx = readsAx or plt.subplot(gs[1, 0])
    else:
        readsAx = readsAx or plt.subplot(111)

    # Make a deep copy of the title alignments. We're potentially going to
    # change the HSP scores, the X axis offsets, etc., and we don't want to
    # interfere with the data we were passed.
    titleAlignments = deepcopy(titlesAlignments[title])

    readsAlignments = titlesAlignments.readsAlignments
    subjectIsNucleotides = readsAlignments.params.subjectIsNucleotides

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
        def adjustOffset(offset):
            return offset

    # It would be more efficient to only walk through all HSPs once and
    # compute these values all at once, but for now this is simple and clear.
    maxY = int(ceil(titleAlignments.bestHsp().score.score))
    minY = int(titleAlignments.worstHsp().score.score)
    maxX = max(hsp.readEndInSubject for hsp in titleAlignments.hsps())
    minX = min(hsp.readStartInSubject for hsp in titleAlignments.hsps())

    if xRange == 'subject':
        # We'll display a graph for the full subject range. Adjust X axis
        # min/max to make sure we cover at least zero to the sequence length.
        maxX = max(titleAlignments.subjectLength, maxX)
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
            '%s (%s)\nLength %d %s, %d read%s, %d HSP%s.' %
            (
                fill(titleAlignments.subjectTitle, 80), accession,
                titleAlignments.subjectLength,
                'nt' if subjectIsNucleotides else 'aa',
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
            createFigure=False, showFigure=False,
            readsAx=ax[row][col], quiet=True, idList=idList, xRange=xRange)

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


def alignmentPanelHTML(titlesAlignments, proteinGenomeDB, sortOn='maxScore',
                       outputDir=None, idList=False, equalizeXAxes=False,
                       xRange='subject', logLinearXAxis=False,
                       logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                       rankScores=False, showFeatures=True,
                       subjectType='protein'):
    """
    Produces an HTML index file in C{outputDir} and a collection of alignment
    graphs and FASTA files to summarize the information in C{titlesAlignments}.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param proteinGenomeDB: A L{dark.civ.proteins.SqliteIndex} instance.
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
    @param subjectType: A C{str} indicating whether the matched subjects are
        'protein' or 'genome'. This is used to determine what accession number
        to extract from the matched title.
    @raise TypeError: If C{outputDir} is C{None}.
    @raise ValueError: If C{outputDir} is None or exists but is not a
        directory or if C{xRange} is not "subject" or "reads".
    """
    if subjectType not in ('protein', 'genome'):
        raise ValueError('subjectType must be either "protein" or "genome".')

    if xRange not in ('subject', 'reads'):
        raise ValueError('xRange must be either "subject" or "reads".')

    if equalizeXAxes:
        raise NotImplementedError('This feature is not yet implemented.')

    titles = titlesAlignments.sortTitles(sortOn)

    if os.access(outputDir, os.F_OK):
        # outputDir exists. Check it's a directory.
        if not S_ISDIR(os.stat(outputDir).st_mode):
            raise ValueError('%r is not a directory.' % outputDir)
    else:
        if outputDir is None:
            raise ValueError('The outputDir needs to be specified.')
        else:
            os.mkdir(outputDir)

    htmlWriter = AlignmentPanelHTMLWriter(outputDir, titlesAlignments)

    getAccession = (proteinGenomeDB.proteinAccession
                    if subjectType == 'protein'
                    else proteinGenomeDB.genomeAccession)

    for title in titles:
        accession = getAccession(title)
        imageBasename = accession + '.png'
        imageFile = join(outputDir, imageBasename)
        graphInfo = alignmentGraph(
            titlesAlignments, title, accession, addQueryLines=True,
            showFeatures=showFeatures, rankScores=rankScores,
            logLinearXAxis=logLinearXAxis, logBase=logBase,
            showFigure=False, imageFile=imageFile,
            quiet=True, idList=idList, xRange=xRange)

        # Close the image plot to make sure memory is flushed.
        plt.close()
        htmlWriter.addImage(imageBasename, accession, title, graphInfo)

    htmlWriter.close()

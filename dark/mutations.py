from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

from dark import ncbidb


def basePlotter(blastHits, title):
    """
    Plot the reads and the subject, so that bases in the reads which are
    different from the subject are shown. Else a '.' is shown.
    like so:
    subject_gi  ATGCGTACGTACGACACC
    read_1         A......TTC..T

    @param blastHits: A L{dark.blast.BlastHits} instance.
    @param title: A C{str} sequence title that was hit by BLAST. We plot the
        reads that hit this title.
    """
    result = []
    params = blastHits.plotParams
    assert params is not None, ('Oops, it looks like you forgot to run '
                                'computePlotInfo.')

    sequence = ncbidb.getSequence(title, blastHits.records.blastDb)
    subject = sequence.seq
    gi = title.split('|')[1]
    sub = '%s \t \t \t %s' % (gi, subject)
    result.append(sub)

    plotInfo = blastHits.titles[title]['plotInfo']
    items = plotInfo['items']
    print 'items', len(items)
    count = 0
    for item in items:
        count += 1
        hsp = item['hsp']
        queryTitle = blastHits.fasta[item['readNum']].id
        # If the product of the subject and query frame values is +ve,
        # then they're either both +ve or both -ve, so we just use the
        # query as is. Otherwise, we need to reverse complement it.
        if item['frame']['subject'] * item['frame']['query'] > 0:
            query = blastHits.fasta[item['readNum']].seq
            reverse = False
        else:
            # One of the subject or query has negative sense.
            #print 'reversing query'
            query = blastHits.fasta[
                item['readNum']].reverse_complement().seq
            reverse = True
        query = query.upper()
        queryStart = hsp['queryStart']
        subjectStart = hsp['subjectStart']
        queryEnd = hsp['queryEnd']
        subjectEnd = hsp['subjectEnd']

        # Before comparing the read to the subject, make a string of the
        # same length as the subject, which contains the read and
        # has ' ' at the offsets.
        # 3 parts need to be taken into account:
        # 1) the left offset (if the query doesn't stick out to the left)
        # 2) the query. if the frame is -1, it has to be reversed.
        # The query consists of 3 parts: left, middle (control for gaps)
        # 3) the right offset

        # do part 1) and 2).
        if queryStart < 0:
            # The query is sticking out to the left
            leftQuery = ''
            if subjectStart == 0:
                # The match starts at the first base of the subject
                middleLeftQuery = ''
                pass
            else:
                # The match starts into the subject.
                # Determine the length of the not matching query
                # part to the left
                leftOffset = -1 * queryStart
                rightOffset = subjectStart + leftOffset
                middleLeftQuery = query[leftOffset:rightOffset]
        else:
            # The query is not sticking out to the left
            # make the left offset
            leftQuery = queryStart * ' '

            leftQueryOffset = subjectStart - queryStart
            middleLeftQuery = query[:leftQueryOffset]

        # do part 3)
        # disregard gaps in subject while adding
        matchQuery = item['origHsp'].query
        matchSubject = item['origHsp'].sbjct
        index = 0
        mid = ''
        for item in range(len(matchQuery)):
            if matchSubject[index] != ' ':
                mid = mid + matchQuery[index]
                index += 1
            else:
                index += 1
        # if the query has been reversed, turn it around
        if reverse:
            rev = ''
            toReverse = mid
            reverseDict = {' ': ' ', '-': '-', 'A': 'T', 'T': 'A',
                           'C': 'G', 'G': 'C', '.': '.', 'N': 'N'}
            for item in toReverse:
                newItem = reverseDict[item]
                rev = rev + newItem
            mid = rev[::-1]

        middleQuery = middleLeftQuery + mid

        # add right not-matching part of the query
        rightQueryOffset = queryEnd - subjectEnd
        rightQuery = query[-rightQueryOffset:]
        middleQuery += rightQuery

        read = leftQuery + middleQuery

        # do part 3)
        sbLen = len(subject)
        offset = sbLen - len(read)
        # if the read is sticking out to the right
        # chop it off
        if offset < 0:
            read = read[:offset]
        # if it's not sticking out, fill the space with ' '
        elif offset > 0:
            addOn = offset * ' '
            read = read + addOn

        # compare the subject and the read, make a string
        # called 'comparison', which contains a '.' if the bases
        # are equal and the letter of the read if they are not.
        index = 0
        comparison = ''
        for base in range(len(read)):
            subjectBase = subject[index]
            readBase = read[index]
            if readBase == ' ':
                comparison += ' '
            elif readBase == subjectBase:
                comparison += '.'
            elif readBase != subjectBase:
                comparison += readBase
            index += 1
        que = '%s \t %s' % (queryTitle, comparison)
        result.append(que)

        # sanity checks
        assert (len(comparison) == len(subject)), (
            '%d "!=" %d' % (len(comparison), len(subject)))

        index = 0
        if comparison[index] == ' ':
            index += 1
        else:
            print index
            start = index - 1
            assert (start == queryStart or start == -1), (start, queryStart)

    return result


def getAPOBECFrequencies(dotAlignment, orig, new, pattern):
    """
    Gets mutation frequencies if they are in a certain pattern.

    @param dotAlignment: result from calling basePlotter
    @param orig: A C{str}, naming the original base
    @param new: A C{str}, what orig was mutated to
    @param pattern: A C{str}m which pattern we're looking for
        (must be one of 'cPattern', 'tPattern')
    """
    cPattern = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT',
                'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT']
    tPattern = ['ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
                'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT']
    # choose the right pattern
    if pattern == 'cPattern':
        patterns = cPattern
        middleBase = 'C'
    else:
        patterns = tPattern
        middleBase = 'T'
    # generate the freqs dict with the right pattern
    freqs = {}
    for pattern in patterns:
        freqs[pattern] = 0
    # get the subject sequence from dotAlignment
    subject = dotAlignment[0].split('\t')[3]
    # exclude the subject from the dotAlignment, so just the queries
    # are left over
    queries = dotAlignment[1:]
    for item in queries:
        query = item.split('\t')[1]
        index = 0
        for queryBase in query:
            qBase = query[index]
            sBase = subject[index]
            if qBase == new and sBase == orig:
                try:
                    plusSb = subject[index + 1]
                    minusSb = subject[index - 1]
                except IndexError:
                    plusSb = 'end'
                motif = '%s%s%s' % (minusSb, middleBase, plusSb)
                if motif in freqs.keys():
                    freqs[motif] += 1
            index += 1

    return freqs


def getCompleteFreqs(blastHits):
    """
    Make a dictionary which collects all mutation frequencies from
    all reads.
    Calls basePlotter to get dotAlignment, which is passed to
    getAPOBECFrequencies with the respective parameter, to collect
    the frequencies.

    @param blastHits: A L{dark.blast.BlastHits} instance.
    """
    allFreqs = {}
    for title in blastHits.titles:
        allFreqs[title] = {
            'C>A': {},
            'C>G': {},
            'C>T': {},
            'T>A': {},
            'T>C': {},
            'T>G': {},
        }
        basesPlotted = basePlotter(blastHits, title)
        for mutation in allFreqs[title]:
            orig = mutation[0]
            new = mutation[2]
            if orig == 'C':
                pattern = 'cPattern'
            else:
                pattern = 'tPattern'
            freqs = getAPOBECFrequencies(basesPlotted, orig, new, pattern)
            allFreqs[title][mutation] = freqs
    return allFreqs


def makeFrequencyGraph(allFreqs, title, substitution, pattern,
                       color='blue', createFigure=True, showFigure=True,
                       readsAx=False):
    """
    For a title, make a graph showing the frequencies.

    @param allFreqs: result from getCompleteFreqs
    @param title: A C{str}, title of virus of which frequencies should be
        plotted.
    @param substitution: A C{str}, which substitution should be plotted;
        must be one of 'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'.
    @param pattern: A C{str}, which pattern we're looking for ( must be
        one of 'cPattern', 'tPattern')
    @param color: A C{str}, color of bars.
    """
    cPattern = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT',
                'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT']
    tPattern = ['ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
                'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT']

    # choose the right pattern
    if pattern == 'cPattern':
        patterns = cPattern
    else:
        patterns = tPattern

    fig = plt.figure(figsize=(10, 10))
    ax = readsAx or fig.add_subplot(111)
    # how many bars
    N = 16
    ind = np.arange(N)
    width = 0.4
    # make a list in the right order, so that it can be plotted easily
    toPlot = allFreqs[title][substitution]
    index = 0
    data = []
    for item in patterns:
        newData = toPlot[patterns[index]]
        data.append(newData)
        index += 1
    # create the bars
    rects1 = ax.bar(ind, data, width, color=color)
    maxY = np.max(data) + 5
    # axes and labels
    if createFigure:
        title = title.split('|')[4][:50]
        ax.set_title('%s \n %s' % (title, substitution), fontsize=20)
        ax.set_ylim(0, maxY)
        ax.set_ylabel('Absolute Number of Mutations', fontsize=16)
        ax.set_xticks(ind+width)
        xTickNames = ax.set_xticklabels(patterns, rotation=45, fontsize=8)
    if not createFigure:
        ax.set_xticks(ind+width)
        xTickNames = ax.set_xticklabels(patterns, rotation=45, fontsize=0)
    if createFigure:
        if showFigure:
            plt.show()
    return maxY


def makeFrequencyPanel(allFreqs, patientName):
    """
    For a title, make a graph showing the frequencies.

    @param allFreqs: result from getCompleteFreqs
    @param patientName: A C{str}, title for the panel
    """
    origMaxY = 0
    cols = 6
    rows = len(allFreqs)
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    substitutions = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    colors = ['blue', 'black', 'red', 'yellow', 'green', 'orange']

    for i, title in enumerate(allFreqs):
        index = 0
        for subst in allFreqs[title]:
            substitution = substitutions[index]
            print i, index, title, 'substitution', substitutions[index]
            if substitution[0] == 'C':
                pattern = 'cPattern'
            else:
                pattern = 'tPattern'
            maxY = makeFrequencyGraph(allFreqs, title, substitution,
                                      pattern, color=colors[index],
                                      createFigure=False, showFigure=False,
                                      readsAx=ax[i][index])
            if maxY > origMaxY:
                origMaxY = maxY

            # add title for individual plot.
            # if used for other viruses, this will have to be adapted.
            if index == 0:
                gi = title.split('|')[1]
                titles = title.split(' ')
                try:
                    typeIndex = titles.index('type')
                except ValueError:
                    typeIndex = -1
                if typeIndex > -1:
                    typeNumber = titles[typeIndex+1]
                else:
                    typeNumber = 'gi: %s' % gi

                ax[i][index].set_ylabel(('Type %s' % (typeNumber)),
                                        fontsize=10)
            # add xAxis tick labels
            if i == 0:
                ax[i][index].set_title(substitution, fontsize=13)
            if i == len(allFreqs)-1 or i == (len(allFreqs)-1)/2:
                if index < 3:
                    pat = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG',
                           'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC',
                           'TCG', 'TCT']
                else:
                    pat = ['ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG',
                           'CTT', 'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC',
                           'TTG', 'TTT']
                ax[i][index].set_xticklabels(pat, rotation=45, fontsize=8)
            index += 1

    # make Y-axis equal
    for i, title in enumerate(allFreqs):
        index = 0
        for substitution in allFreqs[title]:
            a = ax[i][index]
            a.set_ylim([0, origMaxY])
            index += 1
    # add title of whole panel
    figure.suptitle('Mutation Signatures in %s' % patientName, fontsize=20)
    figure.set_size_inches(5 * cols, 3 * rows, forward=True)
    figure.show()

    return allFreqs

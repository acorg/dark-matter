import re
from math import ceil
from IPython.display import HTML
from collections import defaultdict
from random import uniform
from time import ctime, time
import subprocess
from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
from cStringIO import StringIO
from math import log10
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib import gridspec
from urllib2 import URLError

from dark import html
from dark.baseimage import BaseImage
from dark.conversion import readJSONRecords
from dark.dimension import dimensionalIterator
from dark.hsp import normalizeHSP, printHSP
from dark.intervals import OffsetAdjuster, ReadIntervals
from dark import features

Entrez.email = 'tcj25@cam.ac.uk'

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

DEFAULT_LOG_LINEAR_X_AXIS_BASE = 1.1


def readBlastRecords(filename, limit=None):
    """
    Read BLAST records in either XML or JSON format.
    """
    if filename.endswith('.xml'):
        fp = open(filename)
        records = NCBIXML.parse(fp)
    elif filename.endswith('.json'):
        records = readJSONRecords(filename)
        fp = None
    else:
        raise ValueError('Unknown BLAST record file type.')
    for count, record in enumerate(records):
        if limit is not None and count == limit:
            break
        else:
            yield record
    if fp:
        fp.close()


def printBlastRecord(record):
    for key in sorted(record.__dict__.keys()):
        if key not in ['alignments', 'descriptions', 'reference']:
            print '%s: %r' % (key, record.__dict__[key])
    print 'alignments: (%d in total):' % len(record.alignments)
    for i, alignment in enumerate(record.alignments):
        print '  description %d:' % (i + 1)
        for attr in ['accession', 'bits', 'e', 'num_alignments', 'score',
                     'title']:
            print '    %s: %s' % (attr, getattr(record.descriptions[i], attr))
        print '  alignment %d:' % (i + 1)
        for attr in 'accession', 'hit_def', 'hit_id', 'length', 'title':
            print '    %s: %s' % (attr, getattr(alignment, attr))
        print '    HSPs (%d in total):' % len(alignment.hsps)
        for hspIndex, hsp in enumerate(alignment.hsps, start=1):
            print '      hsp %d:' % hspIndex
            printHSP(hsp, '        ')


def summarizeAllRecords(filename, eCutoff=None):
    """
    Read a file of BLAST records and return a dictionary keyed by sequence
    title, with values containing information about the number of times the
    sequence was hit, the e value (from the best HSP), and the sequence
    length.
    eCutoff: A float e value. Hits with e value greater than or equal to
        this will be ignored.
    """
    start = time()
    result = {}
    for record in readBlastRecords(filename):
        for index, alignment in enumerate(record.alignments):
            title = record.descriptions[index].title
            if title in result:
                item = result[title]
            else:
                item = result[title] = {
                    'count': 0,  # Should have been called 'readCount'.
                    'eValues': [],
                    'length': alignment.length,  # This is the subject length.
                    'reads': set(),
                    'title': title,
                }

            if eCutoff is None or alignment.hsps[0].expect <= eCutoff:
                item['count'] += 1
                item['eValues'].append(alignment.hsps[0].expect)
                # record.query is the name of the read in the FASTA file.
                item['reads'].add(record.query)

    # remove subjects with no hits below the cutoff:
    titles = result.keys()
    for title in titles:
        if result[title]['count'] == 0:
            del result[title]
    # compute mean and median evalues:
    for title, value in result.iteritems():
        value['eMean'] = sum(value['eValues']) / float(value['count'])
        value['eMedian'] = np.median(value['eValues'])
        value['reads'] = sorted(value['reads'])

    stop = time()
    report('Record summary generated in %.3f mins.' % ((stop - start) / 60.0))
    return result


def _sortSummary(filenameOrSummary, attr='count', reverse=False):
    """
    Given a filename of BLAST output or a dict that already summarizes
    BLAST output, produce an HTML object with the records sorted by the
    given attribute ('count', 'eMean', 'eMedian', or 'length').
    """
    if isinstance(filenameOrSummary, dict):
        summary = filenameOrSummary
    else:
        summary = summarizeAllRecords(filenameOrSummary)
    if attr not in ('count', 'eMean', 'eMedian', 'length'):
        raise ValueError("attr must be one of 'count', 'eMean', "
                         "'eMedian', or 'length'")
    out = []
    titles = sorted(summary.keys(), key=lambda title: summary[title][attr],
                    reverse=reverse)
    for i, title in enumerate(titles, start=1):
        item = summary[title]
        link = html.NCBISequenceLink(title, title)
        out.append(
            '%3d: count=%4d, len=%7d, median(e)=%20s mean(e)=%20s: %s' %
            (i, item['count'], item['length'], item['eMedian'], item['eMean'],
             link))
    return HTML('<pre><tt>' + '<br/>'.join(out) + '</tt></pre>')


def summarizeAllRecordsByMeanEValueHTML(filenameOrSummary):
    return _sortSummary(filenameOrSummary, 'eMean')


def summarizeAllRecordsByMedianEValueHTML(filenameOrSummary):
    return _sortSummary(filenameOrSummary, 'eMedian')


def summarizeAllRecordsByCountHTML(filenameOrSummary):
    return _sortSummary(filenameOrSummary, 'count', reverse=True)


def summarizeAllRecordsByLengthHTML(filenameOrSummary):
    return _sortSummary(filenameOrSummary, 'length', reverse=True)


def getAllHitsForSummary(summary, recordFilename):
    """
    For the keys of summary (these are sequence titles), pull out a list of
    hit ids and find all those hits in recordFilename.
    """
    titles = sorted(summary.keys())
    hitIds = set([title.split(' ')[0] for title in titles])
    return list(findHits(recordFilename, hitIds))


def filterRecords(summary, filterFunc):
    """
    Pass each of the items in summary (as produced by summarizeAllRecords)
    to a filter function and return a dict subset of those for which the
    filter returns True.
    """
    result = {}
    for title, item in summary.iteritems():
        if filterFunc(item):
            result[title] = item
    return result


def interestingRecords(summary, titleRegex=None, minSequenceLen=None,
                       maxSequenceLen=None, minMatchingReads=None,
                       maxMeanEValue=None, maxMedianEValue=None,
                       negativeTitleRegex=None):
    """
    Given a summary of BLAST results, produced by summarizeAllRecords, return
    a dictionary consisting of just the interesting records.

    summary: the dict output of summarizeAllRecords
    titleRegex: a regex that sequence titles must match.
    negativeTitleRegex: a regex that sequence titles must not match.
    minSequenceLen: sequences of lesser length will be elided.
    maxSequenceLen: sequences of greater length will be elided.
    minMatchingReads: sequences that are matched by fewer reads
        will be elided.
    maxMeanEValue: sequences that are matched with a mean e-value
        that is greater will be elided.
    maxMedianEValue: sequences that are matched with a median e-value
        that is greater will be elided.
    """
    result = {}
    if titleRegex is not None:
        titleRegex = re.compile(titleRegex, re.I)
    if negativeTitleRegex is not None:
        negativeTitleRegex = re.compile(negativeTitleRegex, re.I)
    for title, item in summary.iteritems():
        if titleRegex and titleRegex.search(title) is None:
            continue
        if negativeTitleRegex and negativeTitleRegex.search(title) is not None:
            continue
        if minSequenceLen is not None and item['length'] < minSequenceLen:
            continue
        if maxSequenceLen is not None and item['length'] > maxSequenceLen:
            continue
        if minMatchingReads is not None and item['count'] < minMatchingReads:
            continue
        if maxMeanEValue is not None and item['eMean'] > maxMeanEValue:
            continue
        if maxMedianEValue is not None and item['eMedian'] > maxMedianEValue:
            continue
        result[title] = item
    return result


def getSequence(hitId, db='nt'):
    """
    hitId: the hit_id field from a BLAST record hsp. Of the form
        'gi|63148399|gb|DQ011818.1|' or anything recognized by the -entry param
        of blastdbcmd.
    db: the str name of the BLAST database to search.
    """
    fasta = subprocess.check_output(
        ['blastdbcmd', '-entry', hitId, '-db', db])
    return SeqIO.read(StringIO(fasta), 'fasta')


def findHits(recordFilename, hitIds, limit=None):
    """
    Look in recordFilename (a BLAST XML output file) for hits with ids in the
    set hitIds.

    recordFilename: the str file name to read BLAST records from. Must have
        contents produced via the "-outfmt 5" given on the blast command line.
    hitIds: a set of hit_id field values from a BLAST record hsp. Each hitId is
        of the form 'gi|63148399|gb|DQ011818.1|' or anything recognized by the
        -entry param of blastdbcmd.
    limit: the int number of records to read from recordFilename.

    Return a generator that yields (read number, hit id, hit length, hsps)
    tuples.
    """
    start = time()
    report('Looking for hits on %d sequence ids in %s' %
           (len(hitIds), recordFilename))
    hitCount = 0
    for readNum, record in enumerate(
            readBlastRecords(recordFilename, limit=limit)):
        for alignment in record.alignments:
            if alignment.hit_id in hitIds:
                hitCount += 1
                yield (readNum, alignment.hit_id, alignment.length,
                       alignment.hsps, record.query)
    stop = time()
    report('%d hits found in %.3f mins.' % (hitCount, (stop - start) / 60.0))


def getSeqFromGenbank(hitId):
    """
    hitId: either a hit from a BLAST record, in the form
        'gi|63148399|gb|DQ011818.1|' in which case we use the 2nd field, as
        delimited by '|', to fetch from Genbank.  Or, a gi number (the 2nd
        field just mentioned).

    NOTE: this uses the network!  Also, there is a 3 requests/second limit
    imposed by NCBI on these requests so be careful or your IP will be banned.
    """
    try:
        gi = hitId.split('|')[1]
    except IndexError:
        # Assume we have a gi number directly, and make sure it's a string.
        gi = str(hitId)

    try:
        client = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text',
                               id=gi)
    except URLError:
        return None
    else:
        record = SeqIO.read(client, 'gb')
        client.close()
        return record


def summarizeHits(hits, fastaFilename, eCutoff=None,
                  maxHspsPerHit=None, minStart=None, maxStop=None,
                  logLinearXAxis=False,
                  logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                  randomizeZeroEValues=False):
    """
    Summarize the information found in 'hits'.

    hits: The result of calling findHits (above).
    fastaFilename: The name of the FASTA file containing query sequences,
        or a list of fasta sequences from SeqIO.parse
    eCutoff: A float e value. Hits with e value greater than or equal to
        this will be ignored.
    maxHspsPerHit: The maximum number of HSPs to examine for each hit.
    minStart: Reads that start before this subject offset should not be
        returned.
    maxStop: Reads that end after this subject offset should not be returned.
    logLinearXAxis: if True, convert read offsets so that empty regions in the
        plot we're preparing will only be as wide as their logged actual
        values.
    logBase: The base of the logarithm to use if logLinearXAxis is True.

    Return a 2-tuple: a list of the SeqIO parsed fasta file, and a dict whose
        keys are sequence ids and whose values give hit summary information.
    """

    if isinstance(fastaFilename, str):
        fasta = list(SeqIO.parse(fastaFilename, 'fasta'))
    else:
        fasta = fastaFilename
    zeroEValueUpperRandomIncrement = 150

    def resultDict(sequenceLen):
        result = {
            'hitCount': 0,
            'items': [],
            'maxE': 0.0,
            'minE': 1000,  # Something ridiculously large.
            'maxX': maxStop or sequenceLen,
            'minX': minStart or 0,
            'sequenceLen': sequenceLen,
            'zeroEValueFound': False,
        }

        if logLinearXAxis:
            result['readIntervals'] = ReadIntervals(sequenceLen)
            result['offsetAdjuster'] = None  # Will be set below.

        return result

    result = {}

    # Extract all e values.
    for sequenceId, hitId, hitLen, hsps, query in hits:
        if hitId not in result:
            result[hitId] = resultDict(hitLen)
        hitInfo = result[hitId]
        # Manually count hits. 'hits' may be a generator so we can't use len().
        hitInfo['hitCount'] += 1
        queryLen = len(fasta[sequenceId])

        for hspCount, hsp in enumerate(hsps, start=1):
            if maxHspsPerHit is not None and hspCount > maxHspsPerHit:
                break
            try:
                normalized = normalizeHSP(hsp, queryLen)
            except AssertionError:
                # TODO: Remove these prints, and the surrounding try/except
                # once we're sure we have HSP normalization right.
                #
                # print 'Assertion error in utils calling normalizeHSP'
                # print 'sequenceId: %s' % sequenceId
                # print 'hitId: %s' % hitId
                # print 'hitLen: %s' % hitLen
                # print 'query: %s' % query
                # printHSP(hsp)
                raise
            if ((minStart is not None and normalized['queryStart'] < minStart)
                    or (maxStop is not None and
                        normalized['queryEnd'] > maxStop)):
                continue
            if logLinearXAxis:
                hitInfo['readIntervals'].add(normalized['queryStart'],
                                             normalized['queryEnd'])
            e = hsp.expect
            if eCutoff is not None and e >= eCutoff:
                continue
            if e == 0.0:
                convertedE = None
                hitInfo['zeroEValueFound'] = True
            else:
                convertedE = -1.0 * log10(e)
                if convertedE > hitInfo['maxE']:
                    hitInfo['maxE'] = convertedE
                # Don't use elif for testing minE. Both conditions can be true.
                if convertedE < hitInfo['minE']:
                    hitInfo['minE'] = convertedE
            if normalized['queryStart'] < hitInfo['minX']:
                hitInfo['minX'] = normalized['queryStart']
            if normalized['queryEnd'] > hitInfo['maxX']:
                hitInfo['maxX'] = normalized['queryEnd']
            hitInfo['items'].append({
                'convertedE': convertedE,
                'hsp': normalized,
                'origHsp': hsp,
                'queryLen': queryLen,
                'sequenceId': sequenceId,
                'frame': {
                    'query': hsp.frame[0],
                    'subject': hsp.frame[1],
                },
                'query': query
            })

    for hitInfo in result.itervalues():
        # For each sequence we have hits on, set the expect values that
        # were zero to a randomly high value (higher than the max e value
        # we just calculated).
        maxEIncludingRandoms = hitInfo['maxE']
        if hitInfo['zeroEValueFound']:
            if randomizeZeroEValues:
                for item in hitInfo['items']:
                    if item['convertedE'] is None:
                        item['convertedE'] = e = (
                            hitInfo['maxE'] + 2 + uniform(
                                0, zeroEValueUpperRandomIncrement))
                        if e > maxEIncludingRandoms:
                            maxEIncludingRandoms = e

            else:
                for item in hitInfo['items']:
                    if item['convertedE'] is None:
                        maxEIncludingRandoms += 1
                        item['convertedE'] = maxEIncludingRandoms

        hitInfo['maxEIncludingRandoms'] = maxEIncludingRandoms

        # Adjust all HSPs if we're doing a log/linear X axis.
        if logLinearXAxis:
            adjuster = OffsetAdjuster(hitInfo['readIntervals'], base=logBase)
            minX = maxX = None
            for item in hitInfo['items']:
                adjusted = adjuster.adjustNormalizedHSP(item['hsp'])
                item['hsp'] = adjusted
                if minX is None or adjusted['queryStart'] < minX:
                    minX = adjusted['queryStart']
                if maxX is None or adjusted['queryEnd'] > maxX:
                    maxX = adjusted['queryEnd']

            # Adjust minX and maxX if we have gaps at the start or end of
            # the subject.
            gaps = list(hitInfo['readIntervals'].walk())
            if gaps:
                # Check start of first gap:
                intervalType, (start, stop) = gaps[0]
                if intervalType == ReadIntervals.EMPTY:
                    adjustedStart = adjuster.adjustOffset(start)
                    if adjustedStart < minX:
                        minX = adjustedStart
                # Check stop of last gap:
                intervalType, (start, stop) = gaps[-1]
                if intervalType == ReadIntervals.EMPTY:
                    adjustedStop = adjuster.adjustOffset(stop)
                    if adjustedStop > maxX:
                        maxX = adjustedStop

            hitInfo.update({
                'minX': minX,
                'maxX': maxX,
                'offsetAdjuster': adjuster,
            })

    return fasta, result


def convertSummaryEValuesToRanks(hitInfo):
    """
    Change e values for the reads that hit a sequence to be their ranks.

    hitInfo is a dict of information of BLAST hits against a sequence. It is
    one of the values in the result dictionary returned by summarizeHits.
    """
    result = hitInfo.copy()
    items = result['items']
    items.sort(key=lambda item: item['convertedE'])
    for i, item in enumerate(items, start=1):
        item['convertedE'] = i
    result['maxE'] = result['maxEIncludingRandoms'] = len(items)
    result['minE'] = 1
    result['zeroEValueFound'] = False
    return result


def alignmentGraph(recordFilenameOrHits, hitId, fastaFilename, db='nt',
                   addQueryLines=True, showFeatures=True, eCutoff=1e-2,
                   maxHspsPerHit=None, colorQueryBases=False, minStart=None,
                   maxStop=None, createFigure=True, showFigure=True,
                   readsAx=None, rankEValues=False, imageFile=None,
                   quiet=False, idList=False, xRange='subject',
                   logLinearXAxis=False,
                   logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                   randomizeZeroEValues=False):

    """
    Align a set of BLAST hits against a sequence.

    recordFilenameOrHits: if a string, the BLAST XML output file. Else, a
        dict of BLAST records (e.g., from interestingRecords or
        summarizeAllRecords).
    hitId: the str sequence id to examine the BLAST output for hits against.
    fastaFilename: The name of the FASTA file containing query sequences,
        or a list of fasta sequences from SeqIO.parse
    db: the BLAST db to use to look up the target and, if given, actual
        sequence.
    addQueryLines: if True, draw query lines in full (these will then be partly
        overdrawn by the HSP match against the subject). These are the
        'whiskers' that potentially protrude from each side of a query.
    showFeatures: if True, look online for features of the subject sequence
        (given by hitId).
    eCutoff: e values greater than or equal to this will be ignored.
    maxHspsPerHit: A numeric max number of HSPs to show for each hit on hitId.
    colorQueryBases: if True, color each base of a query string. If True, then
        addQueryLines is meaningless since the whole query is shown colored.
    minStart: Reads that start before this subject offset should not be shown.
    maxStop: Reads that end after this subject offset should not be shown.
    createFigure: If True, create a figure and give it a title.
    showFigure: If True, show the created figure. Set this to False if you're
        creating a panel of figures or just want to save an image (with
        imageFile).
    readsAx: If not None, use this as the subplot for displaying reads.
    rankEValues: If True, display reads with a Y axis coord that is the rank of
        the e value (sorted decreasingly).
    imageFile: If not None, specifies a filename to write the image to.
    quiet: If True, don't print progress / timing output.
    idList: a dictionary. The keys is a color and the values is a list of
        read identifiers that should be colored in the respective color.
    xRange: set to either 'subject' or 'reads' to indicate the range of the
        X axis.
    logLinearXAxis: if True, convert read offsets so that empty regions in the
        plot we're preparing will only be as wide as their logged actual
        values.
    logBase: The base of the logarithm to use if logLinearXAxis is True.
    randomizeZeroEValues: if true places all reads with EValues of 0.0 at
        random positions above the highest read that is not zero. If false,
        places them ranked.
    """

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    start = time()
    sequence = getSequence(hitId, db)

    if createFigure:
        dpi = 80
        # width = max(float(len(sequence)) / float(dpi), 25)
        width = 20
        figure = plt.figure(figsize=(width, 20), dpi=dpi)

    createdReadsAx = readsAx is None

    if showFeatures:
        gbSeq = getSeqFromGenbank(hitId)
        gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 12])
        featureAx = plt.subplot(gs[0, 0])
        orfAx = plt.subplot(gs[1, 0])
        orfReversedAx = plt.subplot(gs[2, 0])
        readsAx = readsAx or plt.subplot(gs[3, 0])
    else:
        featureEndpoints = []
        readsAx = readsAx or plt.subplot(111)

    if isinstance(recordFilenameOrHits, str):
        allhits = findHits(recordFilenameOrHits, set([hitId]))

    else:
        # The recordFilename is actually a dict of hits.
        allhits = recordFilenameOrHits

    fasta, summary = summarizeHits(
        allhits, fastaFilename, eCutoff=eCutoff,
        maxHspsPerHit=maxHspsPerHit, minStart=minStart, maxStop=maxStop,
        logLinearXAxis=logLinearXAxis, logBase=logBase,
        randomizeZeroEValues=randomizeZeroEValues)

    if rankEValues:
        hitInfo = convertSummaryEValuesToRanks(summary[hitId])
    else:
        hitInfo = summary[hitId]

    items = hitInfo['items']
    maxEIncludingRandoms = int(ceil(hitInfo['maxEIncludingRandoms']))
    maxE = int(ceil(hitInfo['maxE']))
    minE = int(hitInfo['minE'])
    maxX = hitInfo['maxX']
    minX = hitInfo['minX']

    # Add light grey vertical rectangles to show the logarithmic gaps. Add
    # these first so that reads will be plotted on top of them. Only draw
    # gaps that are more than 20 pixels wide as we could have millions of
    # tiny gaps for a bacteria and drawing them all will be slow and only
    # serves to make the entire background grey.
    if logLinearXAxis and len(hitInfo['offsetAdjuster'].adjustments()) < 100:
        offsetAdjuster = hitInfo['offsetAdjuster'].adjustOffset
        for (intervalType, interval) in hitInfo['readIntervals'].walk():
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
                query = fasta[item['sequenceId']].seq
            else:
                # One of the subject or query has negative sense.
                query = fasta[item['sequenceId']].reverse_complement().seq
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
            offsetAdjuster = hitInfo['offsetAdjuster'].adjustOffset
        else:
            offsetAdjuster = lambda x: x

        featureEndpoints = features.addFeatures(featureAx, gbSeq, minX, maxX,
                                                offsetAdjuster)
        hitInfo['features'] = featureEndpoints
        if len(featureEndpoints) < 20:
            for fe in featureEndpoints:
                line = Line2D(
                    [fe['start'], fe['start']],
                    [0, maxEIncludingRandoms + 1], color=fe['color'])
                readsAx.add_line(line)
                line = Line2D(
                    [fe['end'], fe['end']],
                    [0, maxEIncludingRandoms + 1], color='#cccccc')
                readsAx.add_line(line)
            features.addORFs(orfAx, sequence.seq, minX, maxX, featureEndpoints,
                             offsetAdjuster)
        else:
            features.addORFs(orfAx, sequence.seq, minX, maxX, [],
                             offsetAdjuster)

        features.addReversedORFs(orfReversedAx,
                                 sequence.reverse_complement().seq,
                                 minX, maxX, offsetAdjuster)

    hitInfo['features'] = featureEndpoints

    # Add the horizontal divider between the highest e value and the randomly
    # higher ones (if any).
    if hitInfo['zeroEValueFound']:
        line = Line2D([minX, maxX], [maxE + 1, maxE + 1], color='#cccccc',
                      linewidth=1)
        readsAx.add_line(line)

    # Titles, axis, etc.
    if createFigure:
        figure.suptitle('%s (length %d, %d hits)' % (
            sequence.description, len(sequence), hitInfo['hitCount']),
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
        for item in hitInfo['items']:
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
        hitInfo['queryMin'] = queryMin
        hitInfo['queryMax'] = queryMax

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
               ((stop - start) / 60.0, len(fasta), len(items)))

    return hitInfo


def alignmentPanel(summary, recordFilenameOrHits, fastaFilename, db='nt',
                   eCutoff=1e-2, maxHspsPerHit=None, minStart=None,
                   maxStop=None, sortOn='eMedian', rankEValues=False,
                   interactive=True, outputDir=None, idList=False,
                   equalizeXAxes=True, xRange='subject',
                   logLinearXAxis=False,
                   logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                   randomizeZeroEValues=False):
    """
    Produces a rectangular panel of graphs that each contain an alignment graph
    against a given sequence.

    summary: the dict output of summarizeAllRecords (or interestingRecords).
        The keys of this dict are the sequence titles that the reads (according
        to the BLAST hits in recordFilename) will be aligned against.
    recordFilenameOrHits: if a string, the BLAST XML output file. Else, a
        dict of BLAST records (e.g., from interestingRecords or
        summarizeAllRecords).
    fastaFilename: The name of the FASTA file containing query sequences,
        or a list of fasta sequences from SeqIO.parse
    db: the BLAST db to use to look up the target and, if given, actual
        sequence.
    eCutoff: e values greater than or equal to this will be ignored.
    maxHspsPerHit: A numeric max number of HSPs to show for each hit on hitId.
    minStart: Reads that start before this subject offset should not be shown.
    maxStop: Reads that end after this subject offset should not be shown.
    sortOn: The attribute to sort subplots on. Either "eMean", "eMedian",
        "title" or "reads"
    rankEValues: If True, display reads with a Y axis coord that is the rank of
        the e value (sorted decreasingly).
    interactive: If True, we are interactive and should display the panel
        using figure.show etc.
    outputDir: If not None, specifies a directory to write an HTML summary to.
    idList: a dictionary. The keys is a color and the values is a list of
        read identifiers that should be colored in the respective color.
    equalizeXAxes: if True, adjust the X axis on each alignment plot to be
        the same.
    xRange: set to either 'subject' or 'reads' to indicate the range of the
        X axis.
    logLinearXAxis: if True, convert read offsets so that empty regions in the
        plot we're preparing will only be as wide as their logged actual
        values.
    logBase: The base of the logarithm to use if logLinearXAxis is True.
    randomizeZeroEValues: if true places all reads with EValues of 0.0 at
        random positions above the highest read that is not zero. If false,
        places them ranked.
    """

    assert xRange in ('subject', 'reads'), (
        'xRange must be either "subject" or "reads".')

    if not (interactive or outputDir):
        raise ValueError('Either interactive or outputDir must be True')

    start = time()
    # Sort titles by mean eValue then title.
    if sortOn == 'eMean':
        titles = sorted(
            summary.iterkeys(),
            key=lambda title: (summary[title]['eMean'], title))
    if sortOn == 'eMedian':
        titles = sorted(
            summary.iterkeys(),
            key=lambda title: (summary[title]['eMedian'], title))
    elif sortOn == 'reads':
        titles = sorted(
            summary.iterkeys(), reverse=True,
            key=lambda title: (len(summary[title]['reads']), title))
    elif sortOn == 'title':
        titles = sorted(summary.iterkeys())
    else:
        raise ValueError('sortOn must be one of "eMean", "eMedian", '
                         '"title" or "reads"')

    if isinstance(fastaFilename, str):
        fasta = list(SeqIO.parse(fastaFilename, 'fasta'))
    else:
        fasta = fastaFilename

    cols = 5
    rows = int(len(titles) / cols) + (0 if len(titles) % cols == 0 else 1)
    figure, ax = plt.subplots(rows, cols, squeeze=False)
    report('Plotting %d titles in %dx%d grid, sorted on %s' %
           (len(titles), rows, cols, sortOn))

    maxEIncludingRandoms = -1
    maxE = -1
    minE = 1000  # Something improbably large.
    maxX = -1
    minX = 1e10  # Something improbably large.
    queryMax = -1
    queryMin = 1e10  # Something improbably large.

    # postProcessInfo will accumulate per-hit information gathered as we
    # make each graph. This will then be used to post-process all the small
    # plots in the panel to make their appearance uniform in various ways.
    postProcessInfo = defaultdict(dict)

    if isinstance(recordFilenameOrHits, str):
        report('Finding all hits in %s.' % recordFilenameOrHits)
        hitIds = set([title.split(' ')[0] for title in titles])
        allhits = list(findHits(recordFilenameOrHits, hitIds))
        report('Found %d hits.' % len(allhits))
    else:
        # recordFilenameOrHits is already the hits we need.
        allhits = recordFilenameOrHits

    if outputDir:
        htmlOutput = html.AlignmentPanelHTML(outputDir, fasta)

    coords = dimensionalIterator((rows, cols))

    for i, title in enumerate(titles):
        row, col = coords.next()
        print '%d: %s %s' % (i, title, html.NCBISequenceLinkURL(title, ''))
        hitId = title.split(' ')[0]
        if interactive:
            hitInfo = alignmentGraph(
                allhits, hitId, fasta, db=db, addQueryLines=True,
                showFeatures=False, eCutoff=eCutoff,
                maxHspsPerHit=maxHspsPerHit, colorQueryBases=False,
                minStart=minStart, maxStop=maxStop, createFigure=False,
                showFigure=False, readsAx=ax[row][col],
                rankEValues=rankEValues, quiet=True, idList=idList,
                xRange=xRange, logLinearXAxis=logLinearXAxis,
                logBase=logBase, randomizeZeroEValues=randomizeZeroEValues)

        if outputDir:
            imageBasename = '%d.png' % i
            imageFile = '%s/%s' % (outputDir, imageBasename)
            hitInfo = alignmentGraph(
                allhits, hitId, fasta, db=db, addQueryLines=True,
                showFeatures=False, eCutoff=eCutoff,
                maxHspsPerHit=maxHspsPerHit, colorQueryBases=False,
                minStart=minStart, maxStop=maxStop, showFigure=False,
                rankEValues=rankEValues, imageFile=imageFile, quiet=True,
                idList=idList, xRange=xRange, logLinearXAxis=logLinearXAxis,
                logBase=logBase, randomizeZeroEValues=randomizeZeroEValues)

            # Close the image plot, otherwise it will be displayed when we
            # call plt.show below.
            plt.close()
            htmlOutput.addImage(imageBasename, title, hitInfo)

        # Remember info required for post processing of the whole panel
        # once we've made all the alignment graphs.
        post = postProcessInfo[(row, col)]
        if hitInfo['zeroEValueFound']:
            post['maxE'] = hitInfo['maxE']
        post['maxEIncludingRandoms'] = hitInfo['maxEIncludingRandoms']
        post['maxX'] = hitInfo['maxX']
        post['sequenceLen'] = hitInfo['sequenceLen']
        if logLinearXAxis:
            post['offsetAdjuster'] = hitInfo['offsetAdjuster']
            post['readIntervals'] = hitInfo['readIntervals']

        #calculate eMedian and eMean:
        eValues = [item['convertedE'] for item in hitInfo['items']]

        if rankEValues:
            meanE = medianE = 'ranked'
        elif len(eValues) == 0:
            meanE = medianE = 'none'
        else:
            numberMeanE = sum(eValues) / float(len(eValues))
            numberMedianE = np.median(eValues)
            meanE = '1e-%d' % numberMeanE
            medianE = '1e-%d' % numberMedianE
        ax[row][col].set_title(
            '%d: %s\n%d HSPs shown (%d reads in total)\n%s median, %s mean' % (
                i, title.split(' ', 1)[1][:40], len(hitInfo['items']),
                len(summary[title]['reads']), medianE, meanE), fontsize=10)

        if hitInfo['maxEIncludingRandoms'] > maxEIncludingRandoms:
            maxEIncludingRandoms = hitInfo['maxEIncludingRandoms']
        if hitInfo['maxE'] > maxE:
            maxE = hitInfo['maxE']
        if hitInfo['minE'] < minE:
            minE = hitInfo['minE']
        if hitInfo['maxX'] > maxX:
            maxX = hitInfo['maxX']
        if hitInfo['minX'] < minX:
            minX = hitInfo['minX']
        if xRange == 'reads':
            if hitInfo['queryMin'] < queryMin:
                queryMin = hitInfo['queryMin']
            if hitInfo['queryMax'] > queryMax:
                queryMax = hitInfo['queryMax']

    # Post-process graphs to adjust axes, etc.
    coords = dimensionalIterator((rows, cols))
    for row, col in coords:
        a = ax[row][col]
        a.set_ylim([0, maxEIncludingRandoms + 1])
        if equalizeXAxes:
            if xRange == 'subject':
                a.set_xlim([minX, maxX])
            else:
                a.set_xlim([queryMin, queryMax])
        a.set_yticks([])
        a.set_xticks([])
        hitInfo = postProcessInfo[(row, col)]
        if hitInfo:
            # The graph is non-blank (the last graphs in a panel can be blank).
            if equalizeXAxes and 'maxE' in hitInfo:
                # Overdraw the horizontal divider between the highest e value
                # and the randomly higher ones (if any). We need to do this
                # as the plots will be changing width, to all be as wide as
                # the widest.
                e = hitInfo['maxE']
                line = Line2D([minX, maxX], [e + 1, e + 1], color='#cccccc',
                              linewidth=1)
                a.add_line(line)

            # Add a vertical line at x=0 so we can see the 'whiskers' of
            # reads that extend to the left of the sequence we're aligning
            # against.
            line = Line2D([0, 0], [0, maxEIncludingRandoms + 1],
                          color='#cccccc', linewidth=1)
            a.add_line(line)
            # Add a line on the right of each sub-plot so we can see where
            # the sequence ends (as all panel graphs have the same width and
            # we otherwise couldn't tell).
            sequenceLen = hitInfo['sequenceLen']
            if logLinearXAxis:
                offsetAdjuster = hitInfo['offsetAdjuster'].adjustOffset
                sequenceLen = offsetAdjuster(sequenceLen)
            line = Line2D([sequenceLen, sequenceLen],
                          [0, maxEIncludingRandoms + 1],
                          color='#cccccc', linewidth=1)
            a.add_line(line)

            # Add light grey vertical rectangles to show the logarithmic
            # gaps. Add these only in the region above the highest read in
            # the individual plot. If we simply added the bar to the full
            # height of the plot it would obscure the reads it overlapped.
            if (logLinearXAxis and
                    len(hitInfo['offsetAdjuster'].adjustments()) < 100):
                thisMaxEIncludingRandoms = hitInfo['maxEIncludingRandoms']
                for (intervalType, interval) in hitInfo[
                        'readIntervals'].walk():
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


def report(msg):
    print '%s: %s' % (ctime(time()), msg)


def evalueGraph(records, rows, cols, find=None, titles=True, minHits=1,
                width=5, height=5):
    """
    Produces a rectangular panel of graphs that each show sorted e-values for
    a read. Read hits against a certain strain (see find, below) are
    highlighted.

    find: A function that can be passed a sequence title. If the function
    returns True a (currently) red dot is put into the graph at that point.
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

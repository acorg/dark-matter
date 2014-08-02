from math import log10
import numpy as np
from random import uniform

from dark.blast.hsp import normalizeHSP
from dark.intervals import OffsetAdjuster, ReadIntervals
from dark import mysql

DEFAULT_LOG_LINEAR_X_AXIS_BASE = 1.1


class BlastHits(object):
    """
    Maintain information about a set of sequences hit by a BLAST run.

    @param records: A L{BlastRecords} instance.
    @param readSetFilter: A L{dark.filter.ReadSetFilter} instance or C{None}.
        If not C{None}, this is the read set filter that was used to filter
        the hits for this read set.
    """

    def __init__(self, records, readSetFilter=None):
        self.records = records
        self.readSetFilter = readSetFilter
        self.titles = {}
        self.fasta = None  # Computed (once) in summarizeHits.
        self.plotParams = None  # Set in computePlotInfo.

    def _convertBitScoresToRanks(self, plotInfo):
        """
        Change the bit scores for the reads that hit each sequence to be their
        ranks.

        @param plotInfo: A C{dict} of plot information, as returned by
            C{plotInfoDict} in C{computePlotInfo}.
        """
        items = plotInfo['items']
        items.sort(key=lambda item: item['bitScore'])
        for rank, item in enumerate(items, start=1):
            item['bitScore'] = rank
        plotInfo['bitScoreMax'] = len(items)
        plotInfo['bitScoreMin'] = 1

    def computePlotInfo(self,
                        logLinearXAxis=False,
                        logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                        randomizeZeroEValues=True, rankValues=False):
        """
        Read detailed HSP information about the hit titles in C{self.titles}
        and compute summary statistics on it. The various parameters allow
        us to restrict and transform the data that is read.

        @param logLinearXAxis: if True, convert read offsets so that empty
            regions in the plot we're preparing will only be as wide as their
            logged actual values.
        @param logBase: The base of the logarithm to use if logLinearXAxis is
            C{True}.
        @param randomizeZeroEValues: If C{True}, e-values that are zero will
            be set to a random (extremely good) value.
        @param: rankValues: If C{True}, change the e-values and bit scores for
            the reads for each title to be their rank (worst to best).
        """

        # Reset the plot info for each hit title because we may be called
        # multiple times, with different parameters.
        for title, hitInfo in self.titles.iteritems():
            hitInfo['plotInfo'] = None

        # Save the plotting parameters in use so that things we pass self
        # to, such as alignmentPlot and alignmentPanel, can discover how
        # the data was filtered and summarized etc.
        self.plotParams = {
            'logBase': logBase,
            'logLinearXAxis': logLinearXAxis,
            'randomizeZeroEValues': randomizeZeroEValues,
            'rankValues': rankValues,
        }

        # Read in the read FASTA file if we haven't already.
        if self.fasta is None:
            self.fasta = self._readFASTA()

        zeroEValueUpperRandomIncrement = 150

        def plotInfoDict(sequenceLen):
            """
            Retun a new C{dict} suitable for holding plot information
            for a hit sequence title.

            @param sequenceLen: The C{int} base pairs length of the sequence.
            @return: A C{dict} to hold plot info for a title.
            """
            result = {
                'bitScores': [],
                'hspTotal': 0,  # The total number of HSPs for this title.
                'items': [],  # len() of this = the # of HSPs kept for display.
                'maxE': None,
                'minE': None,
                'maxX': maxStop or sequenceLen,
                'minX': minStart or 0,
                'zeroEValueFound': False,
                'originalEValues': [],
                'readNums': set(),
                # these two below should go away
                'eMean': None,
                'eMedian': None,
            }

            if logLinearXAxis:
                result['readIntervals'] = ReadIntervals(sequenceLen)
                result['offsetAdjuster'] = None  # Will be set below.

            return result

        # The command-line program that was used to run Blast. Will be set
        # below once we have started to read the BLAST records.
        blastApplication = None

        for title, readNum, hsps in self._getHsps():
            if blastApplication is None:
                blastApplication = self.records.blastParams[
                    'application'].lower()
            if self.titles[title]['plotInfo'] is None:
                sequenceLen = self.titles[title]['length']
                self.titles[title]['plotInfo'] = plotInfoDict(sequenceLen)
            plotInfo = self.titles[title]['plotInfo']
            queryLen = len(self.fasta[readNum])
            plotInfo['hspTotal'] += len(hsps)

            for hspCount, hsp in enumerate(hsps, start=1):
                try:
                    normalized = normalizeHSP(hsp, queryLen, blastApplication)
                except AssertionError:
                    # TODO: Remove these prints, and the surrounding try/except
                    # once we're sure we have HSP normalization right.
                    #
                    # print 'Assertion error in utils calling normalizeHSP'
                    # print 'readNum: %s' % readNum
                    # print 'subjectLength: %s' % self.titles[title]['length']
                    # print 'query: %s' % query
                    # from dark.hsp import printHSP
                    # printHSP(hsp)
                    raise

                if logLinearXAxis:
                    plotInfo['readIntervals'].add(normalized['queryStart'],
                                                  normalized['queryEnd'])

                # We are committed to adding this item to plotInfo['items'].
                # Don't add any 'continue' statements below this point.

                e = hsp.expect

                # Retain original evalue below cutoff for eMean and eMedian
                # calculation.
                plotInfo['originalEValues'].append(e)

                # Keep track of the read numbers that are still valid for
                # the plot. This will be a subset of the set of all read
                # numbers that matched the hit; here we're collecting the
                # read numbers that are still valid once we've done further
                # filtering (as just above).
                plotInfo['readNums'].add(readNum)

                if e == 0.0:
                    convertedE = None
                    plotInfo['zeroEValueFound'] = True
                else:
                    convertedE = -1.0 * log10(e)
                    if plotInfo['minE'] is None:
                        # This is the first item being added. Set minE and
                        # maxE unconditionally.
                        # TODO: Delete the following sanity check.
                        assert plotInfo['maxE'] is None
                        plotInfo['minE'] = plotInfo['maxE'] = convertedE
                    else:
                        if convertedE < plotInfo['minE']:
                            plotInfo['minE'] = convertedE
                        elif convertedE > plotInfo['maxE']:
                            plotInfo['maxE'] = convertedE
                plotInfo['bitScores'].append(hsp.bits)
                if normalized['queryStart'] < plotInfo['minX']:
                    plotInfo['minX'] = normalized['queryStart']
                if normalized['queryEnd'] > plotInfo['maxX']:
                    plotInfo['maxX'] = normalized['queryEnd']

                plotInfo['items'].append({
                    'bitScore': hsp.bits,
                    'convertedE': convertedE,
                    'hsp': normalized,
                    'origHsp': hsp,
                    'queryLen': queryLen,
                    'readNum': readNum,
                    'frame': {
                        'query': hsp.frame[0],
                        'subject': hsp.frame[1],
                    },
                })

        for title in self.titles.iterkeys():
            plotInfo = self.titles[title]['plotInfo']
            if plotInfo is None:
                continue
            if not plotInfo['items']:
                # No items were added above for this title. Reset its plotInfo.
                self.titles[title]['plotInfo'] = None
                continue

            # TODO: do something like the following for the e-values.
            bitScores = plotInfo['bitScores']
            plotInfo['bitScoreMin'] = np.min(bitScores)
            plotInfo['bitScoreMax'] = np.max(bitScores)
            plotInfo['bitScoreMean'] = np.mean(bitScores)
            plotInfo['bitScoreMedian'] = np.median(bitScores)
            del plotInfo['bitScores']

            # If plotInfo['minE'] is None, plotInfo['maxE'] will be too. This
            # indicates that all the qualifying e-values found above were zero.
            # We must have found some values (because plotinfo is not None).
            # We can safely set minE and maxE to harmless (zero) values here.
            # Because a zero e-value was found, maxEIncludingRandoms will be
            # set to a higher-than-zero value below and things will work.
            if plotInfo['minE'] is None:
                # TODO: Remove the asserts once we're sure the logic is right.
                # A couple of sanity checks, for now.
                assert plotInfo['maxE'] is None, (
                    "MaxE is %r" % (plotInfo['maxE'],))
                assert plotInfo['zeroEValueFound']
                plotInfo['minE'] = plotInfo['maxE'] = 0

            if rankValues:
                self._convertEValuesToRanks(plotInfo)
                self._convertBitScoresToRanks(plotInfo)

            # For each sequence we have hits on, set the expect values that
            # were zero. If randomizeZeroEValues is True, we set them to
            # randomly high values (higher than the max e value we just
            # calculated). Otherwise, we give them incrementally higher
            # values (like ranking them) above the maximum e value found
            # above.
            maxEIncludingRandoms = plotInfo['maxE']
            if plotInfo['zeroEValueFound']:
                if randomizeZeroEValues:
                    for item in plotInfo['items']:
                        if item['convertedE'] is None:
                            item['convertedE'] = e = (
                                plotInfo['maxE'] + 2 + uniform(
                                    0, zeroEValueUpperRandomIncrement))
                            if e > maxEIncludingRandoms:
                                maxEIncludingRandoms = e
                else:
                    for item in plotInfo['items']:
                        if item['convertedE'] is None:
                            maxEIncludingRandoms += 1
                            item['convertedE'] = maxEIncludingRandoms

            plotInfo['maxEIncludingRandoms'] = maxEIncludingRandoms

            # Adjust all HSPs if we're doing a log/linear X axis.
            if logLinearXAxis:
                adjuster = OffsetAdjuster(plotInfo['readIntervals'],
                                          base=logBase)
                minX = maxX = None
                for item in plotInfo['items']:
                    adjusted = adjuster.adjustNormalizedHSP(item['hsp'])
                    item['hsp'] = adjusted
                    if minX is None or adjusted['queryStart'] < minX:
                        minX = adjusted['queryStart']
                    if maxX is None or adjusted['queryEnd'] > maxX:
                        maxX = adjusted['queryEnd']

                # Adjust minX and maxX if we have gaps at the start or end of
                # the subject.
                gaps = list(plotInfo['readIntervals'].walk())
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

                # Update the plot info with the new min & max X values. We need
                # to check that we actually have new values, as
                # plotInfo['items'] might have been empty.
                if minX is not None:
                    assert maxX is not None  # Sanity check.
                    plotInfo['minX'] = minX
                    plotInfo['maxX'] = maxX

                plotInfo['offsetAdjuster'] = adjuster

            # Calculate eMedian and eMean
            originalEValues = plotInfo['originalEValues']
            plotInfo['originalEMin'] = np.min(originalEValues)
            plotInfo['originalEMean'] = np.mean(originalEValues)
            plotInfo['originalEMedian'] = np.median(originalEValues)
            del plotInfo['originalEValues']

    def getTaxIDFromMySql(self, db=None):
        """
        For each title in C{self.titles}, read the corresponding taxId from
        the gi_taxid_nucl table in a MySQL database called ncbi_taxonomy.
        """
        if db is None:
            openedDb = True
            db = mysql.getDatabaseConnection()
        else:
            openedDb = False
        cursor = db.cursor()

        # for each title (=gi number) get the taxId from the database
        # and add it to self.titles.
        for title in self.titles:
            giNr = int(title.split('|')[1])
            question = 'SELECT taxID from gi_taxid_nucl where gi = %d' % giNr
            cursor.execute(question)
            try:
                result = cursor.fetchone()[0]
            except TypeError:
                result = 'No taxID found'
            self.titles[title]['taxID'] = result

        cursor.close()
        if openedDb:
            db.close()

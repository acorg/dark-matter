from math import log10
import re
import numpy as np
from random import uniform
from IPython.display import HTML
from Bio.Blast import NCBIXML
from Bio import SeqIO

from dark.conversion import readJSONRecords
from dark.hsp import printHSP, normalizeHSP
from dark.html import NCBISequenceLink
from dark.intervals import OffsetAdjuster, ReadIntervals
from dark.simplify import simplifyTitle

DEFAULT_LOG_LINEAR_X_AXIS_BASE = 1.1


def printBlastRecord(record):
    """
    Print a BLAST record.

    @param record: A BioPython C{Bio.Blast.Record.Blast} instance.
    """
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


class BlastRecords(object):
    """
    Hold information about a set of BLAST records.

    @param blastFilename: the C{str} file name containing the BLAST output.
        This can either be an XML (-outfmt 5) BLAST output file or our
        smaller converted JSON equivalent (produced by
        bin/convert-blast-xml-to-json.py from a BLAST XML file).
    @param blastDb: the BLAST database used.
    @param fastaFilename: the C{str} file name containing the sequences that
        were given to BLAST as queries.
    @param limit: An C{int} limit on the number of records to read.
    """

    def __init__(self, blastFilename, fastaFilename, blastDb, limit=None):
        self.blastFilename = blastFilename
        self.fastaFilename = fastaFilename
        self.blastDb = blastDb
        self.limit = limit
        self._length = 0

    def records(self):
        """
        Extract all BLAST records (up to C{self.limit}, if not C{None}).

        @return: A generator that yields BioPython C{Bio.Blast.Record.Blast}
            instances.
        """
        if self.blastFilename.endswith('.xml'):
            fp = open(self.blastFilename)
            records = NCBIXML.parse(fp)
        elif self.blastFilename.endswith('.json'):
            records = readJSONRecords(self.blastFilename)
            fp = None
        else:
            raise ValueError('Unknown BLAST record file type.')

        # Read the records, observing any given limit. There's a little code
        # duplication here in the name of speed - we don't want to repeatedly
        # test if self.limit is None in a loop.
        if self.limit is None:
            for count, record in enumerate(records, start=1):
                yield record
        else:
            limit = self.limit
            for count, record in enumerate(records, start=1):
                if count > limit:
                    count = limit
                    break
                else:
                    yield record

        try:
            self._length = count
        except NameError:
            # If there were no records, count wont even be defined.
            self._length = 0

        if fp:
            fp.close()

    def __len__(self):
        """
        Return the number of BLAST records.

        NOTE: this will return zero if called before C{records()} has been
              called to read the records.

        @return: an C{int}, the number of hits in this set.
        """
        return self._length

    def hits(self):
        """
        Read the BLAST records and return a L{BlastHits} instance.

        @return: A L{BlastHits} instance.
        """
        result = {}
        for readNum, record in enumerate(self.records()):
            for index, alignment in enumerate(record.alignments):
                title = record.descriptions[index].title
                if title in result:
                    hitInfo = result[title]
                else:
                    hitInfo = result[title] = {
                        'readCount': 0,
                        'eValues': [],
                        'length': alignment.length,  # The hit sequence length.
                        'readNums': set(),
                        # 'title': title,  # Not needed, I don't think.
                    }

                hitInfo['readCount'] += 1
                hitInfo['eValues'].append(alignment.hsps[0].expect)
                hitInfo['readNums'].add(readNum)

        blastHits = BlastHits(self)

        # Compute summary stats on e-values and add the hit info to our
        # final result.
        for title, hitInfo in result.iteritems():
            eValues = hitInfo['eValues']
            hitInfo['eMean'] = sum(eValues) / float(hitInfo['readCount'])
            hitInfo['eMedian'] = np.median(eValues)
            hitInfo['eMin'] = min(eValues)
            # Remove the e-values now that we've summarized them.
            del hitInfo['eValues']

            blastHits.addHit(title, hitInfo)

        return blastHits


class BlastHits(object):
    """
    Maintain information about a set of sequences hit by a BLAST run.

    @param records: A L{BlastRecords} instance.
    """

    def __init__(self, records):
        self.records = records
        self.titles = {}
        self.fasta = None  # Computed (once) in summarizeHits.
        self.plotParams = None  # Set in computePlotInfo.

    def __len__(self):
        """
        Return the number of hits.

        @return: an C{int}, the number of hits in this set.
        """
        return len(self.titles)

    def addHit(self, title, hitInfo):
        """
        Add information about a hit against a sequence with title C{title}. If
        the title has already been seen, raise C{KeyError}.

        @param title: A C{str} sequence title.
        @param hitInfo: A C{dict} of information associated with hit against
            the sequence with this title. The dict must contain:

            readCount: the number of reads that hit this sequence.
            length: the length of the sequence.
            reads: a sorted C{list} of read names (as found in e.g., a FASTA
                file).
            eMean: the mean of the e-values of all hits against this sequence.
            eMedian: the median of the e-values of all hits against this
                sequence.
            eMin: the minimum of the e-values of all hits against this
                sequence.
        """

        if title in self.titles:
            raise KeyError('Title %r already present' % (title,))

        self.titles[title] = hitInfo

    def sortTitles(self, by):
        """
        Sort titles by a given attribute and then by title.

        @param by: A C{str}, one of 'eMean', 'eMedian', 'readCount', 'title'.
        @return: A sorted C{list} of titles.
        """
        if by == 'eMean':
            titles = sorted(
                self.titles.iterkeys(),
                key=lambda title: (self.titles[title]['eMean'], title))
        elif by == 'eMedian':
            titles = sorted(
                self.titles.iterkeys(),
                key=lambda title: (self.titles[title]['eMedian'], title))
        elif by == 'readCount':
            titles = sorted(
                self.titles.iterkeys(), reverse=True,
                key=lambda title: (self.titles[title]['readCount'], title))
        elif by == 'length':
            titles = sorted(
                self.titles.iterkeys(), reverse=True,
                key=lambda title: (self.titles[title]['length'], title))
        elif by == 'title':
            titles = sorted(self.titles.iterkeys())
        else:
            raise ValueError('sort attribute must be one of "eMean", '
                             '"eMedian", "title" or "reads".')
        return titles

    def _sortHTML(self, by):
        """
        Return an HTML object with the hits sorted by the given attribute.

        @param by: A C{str}, one of 'eMean', 'eMedian', 'readCount', 'title'.
        @return: An HTML instance with sorted titles and information about
            hit read count, length, and e-values.
        """
        out = []
        for i, title in enumerate(self.sortTitles(by), start=1):
            hitInfo = self.titles[title]
            link = NCBISequenceLink(title, title)
            out.append(
                '%3d: count=%4d, len=%7d, median(e)=%20s mean(e)=%20s: %s' %
                (i, hitInfo['readCount'], hitInfo['length'],
                 hitInfo['eMedian'], hitInfo['eMean'], link))
        return HTML('<pre><tt>' + '<br/>'.join(out) + '</tt></pre>')

    def summarizeByMeanEValueHTML(self):
        return self._sortHTML('eMean')

    def summarizeByMedianEValueHTML(self):
        return self._sortHTML('eMedian')

    def summarizeByCountHTML(self):
        return self._sortHTML('readCount', reverse=True)

    def summarizeByLengthHTML(self):
        return self._sortHTML('length', reverse=True)

    def interesting(self, titleRegex=None, minSequenceLen=None,
                    maxSequenceLen=None, minMatchingReads=None,
                    maxMeanEValue=None, maxMedianEValue=None,
                    negativeTitleRegex=None, maxMinEValue=None,
                    truncateTitlesAfter=None):
        """
        Produce a new L{BlastHits} instance consisting of just the interesting
        hits, as given by our parameters.

        IMPORTANT: the hitInfo entries in the new instance are the same objects
        we have in self.titles. So if our caller changes the result we return,
        they will be changing self too. This could be prevented by making a
        copy of each hitInfo element we have before passing it to result.addHit
        below.  Sharing the hitInfo is usually desirable because it's summary
        information about a BLAST run and so should not be changing, and it
        saves memory.

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
        maxMinEValue: if the minimum e-value for a hit is higher than this
            value, elide the hit. E.g., suppose we are passed a value of
            1e-20, then we should reject any hit whose minimal (i.e., best)
            e-value is bigger than 1e-20. So a hit with minimal e-value of
            1e-10 would not be reported, whereas a hit with a minimal e-value
            of 1e-30 would be.
        truncateTitlesAfter: specify a string that titles will be truncated
            beyond. If a truncated title has already been seen, that title will
            be elided.
        """
        result = BlastHits(self.records)

        if truncateTitlesAfter:
            truncatedTitles = set()
        if titleRegex is not None:
            titleRegex = re.compile(titleRegex, re.I)
        if negativeTitleRegex is not None:
            negativeTitleRegex = re.compile(negativeTitleRegex, re.I)

        for title, hitInfo in self.titles.iteritems():
            if (minSequenceLen is not None and
                    hitInfo['length'] < minSequenceLen):
                continue
            if (maxSequenceLen is not None and
                    hitInfo['length'] > maxSequenceLen):
                continue
            if (minMatchingReads is not None and
                    hitInfo['readCount'] < minMatchingReads):
                continue
            if maxMeanEValue is not None and hitInfo['eMean'] > maxMeanEValue:
                continue
            if (maxMedianEValue is not None and
                    hitInfo['eMedian'] > maxMedianEValue):
                continue
            if maxMinEValue is not None and hitInfo['eMin'] > maxMinEValue:
                continue
            if truncateTitlesAfter:
                # Titles start with something like gi|525472786|emb|HG313807.1|
                # that we need to skip.
                titleSansId = title.split(' ', 1)[1]
                truncated = simplifyTitle(titleSansId, truncateTitlesAfter)
                if truncated in truncatedTitles:
                    # We've already seen this (truncated) title. Skip it.
                    continue
                truncatedTitles.add(truncated)
            # Do the title regex tests last, since they are slowest.
            if titleRegex and titleRegex.search(title) is None:
                continue
            if (negativeTitleRegex and
                    negativeTitleRegex.search(title) is not None):
                continue

            result.addHit(title, hitInfo)

        return result

    def _getHsps(self):
        """
        Extract detailed HSP information for our titles from the BLAST records.

        @return: A generator yielding tuples of the form
            (title, readNumber, hsps)
            Where:
                title is the title of the sequence that was hit.
                readNumber is the index of the read in the records (and hence
                    in the original FASTA file).
                hsps is the list of HSPs for the read matching against the hit.
        """
        titles = self.titles
        for readNum, record in enumerate(self.records.records()):
            for index, alignment in enumerate(record.alignments):
                title = record.descriptions[index].title
                if title in titles:
                    yield (title, readNum, alignment.hsps)

    def _convertEValuesToRanks(self, plotInfo):
        """
        Change e-values for the reads that hit each sequence to be their ranks.

        @param plotInfo: A C{dict} of plot information, as returned by
            C{plotInfoDict} in C{computePlotInfo} and built up there.
        """
        items = plotInfo['items']
        items.sort(key=lambda item: item['convertedE'])
        for rank, item in enumerate(items, start=1):
            item['convertedE'] = rank
        plotInfo['maxE'] = plotInfo['maxEIncludingRandoms'] = len(items)
        plotInfo['minE'] = 1
        plotInfo['zeroEValueFound'] = False

    def _readFASTA(self):
        """
        Read the FASTA data and check its length is compatible with the
        number of BLAST records.

        @return: A C{list} of BioPython sequences.
        """
        fasta = list(SeqIO.parse(self.records.fastaFilename, 'fasta'))
        # Sanity check that the number of reads in the FASTA file is
        # compatible with the number of records in the BLAST file.
        if self.records.limit is None:
            assert len(fasta) == len(self.records), (
                'Sanity check failed: mismatched BLAST and FASTA files. '
                'BLAST file %r contains %d records, whereas FASTA file '
                '%r contains %d sequences.' %
                (self.records.blastFilename, len(self.records),
                 self.records.fastaFilename, len(fasta)))
        else:
            assert len(fasta) >= len(self.records), (
                'Sanity check failed: mismatched BLAST and FASTA files. '
                'BLAST file %r contains at least %d records, whereas '
                'FASTA file %r only contains %d sequences.' %
                (self.records.blastFilename, self.records.limit,
                 self.records.fastaFilename, len(fasta)))
            # Truncate the FASTA to match the limit we have on the
            # number of BLAST records.
            fasta = fasta[:self.records.limit]
        return fasta

    def computePlotInfo(self, eCutoff=None, maxHspsPerHit=None,
                        minStart=None, maxStop=None, logLinearXAxis=False,
                        logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                        randomizeZeroEValues=False, rankEValues=False):
        """
        Summarize the information found in 'hits'.

        @param eCutoff: A float e value. Hits with e value greater than or
            equal to this will be ignored.
        @param maxHspsPerHit: The maximum number of HSPs to examine for each
            hit.
        @param minStart: Reads that start before this subject offset should
            not be returned.
        @param maxStop: Reads that end after this subject offset should not
            be returned.
        @param logLinearXAxis: if True, convert read offsets so that empty
            regions in the plot we're preparing will only be as wide as their
            logged actual values.
        @param logBase: The base of the logarithm to use if logLinearXAxis is
            C{True}.
        @param randomizeZeroEValues: If C{True}, e-values that are zero will
            be set to a random (extremely good) value.
        @param: rankEValues: If C{True}, change the e-values for the reads
            for each title to be their rank (sorted decreasingly).

        @return: A L{BlastHitsSummary} instance.
        """

        # Reset the plot info for each hit title because we may be called
        # multiple times, with different parameters.
        for title, hitInfo in self.titles.iteritems():
            hitInfo['plotInfo'] = None

        # Save the plotting parameters in use so that things we pass self
        # to such as alignmentPlot or alignmentPanel can discover how the
        # data has been filtered and summarized etc.
        self.plotParams = {
            'eCutoff': eCutoff,
            'logBase': logBase,
            'logLinearXAxis': logLinearXAxis,
            'maxHspsPerHit': maxHspsPerHit,
            'maxStop': maxStop,
            'minStart': minStart,
            'randomizeZeroEValues': randomizeZeroEValues,
            'rankEValues': rankEValues,
        }

        # Read in the read FASTA file if we haven't already.
        if self.fasta is None:
            self.fasta = self._readFASTA()

        zeroEValueUpperRandomIncrement = 150

        def plotInfoDict(sequenceLen):
            result = {
                'hspTotal': 0,  # The total number of HSPs for this title.
                'items': [],  # len() of this = the # of HSPs kept for display.
                'maxE': 0.0,
                'minE': 1000,  # Something ridiculously large.
                'maxX': maxStop or sequenceLen,
                'minX': minStart or 0,
                'zeroEValueFound': False,
            }

            if logLinearXAxis:
                result['readIntervals'] = ReadIntervals(sequenceLen)
                result['offsetAdjuster'] = None  # Will be set below.

            return result

        for title, readNum, hsps in self._getHsps():
            if self.titles[title]['plotInfo'] is None:
                sequenceLen = self.titles[title]['length']
                self.titles[title]['plotInfo'] = plotInfoDict(sequenceLen)
            plotInfo = self.titles[title]['plotInfo']
            queryLen = len(self.fasta[readNum])
            plotInfo['hspTotal'] += len(hsps)

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
                    # print 'readNum: %s' % readNum
                    # print 'hitLen: %s' % self.titles[title]['length']
                    # print 'query: %s' % query
                    # from dark.hsp import printHSP
                    # printHSP(hsp)
                    raise
                if ((minStart is not None and
                     normalized['queryStart'] < minStart)
                    or (maxStop is not None and
                        normalized['queryEnd'] > maxStop)):
                    continue
                if logLinearXAxis:
                    plotInfo['readIntervals'].add(normalized['queryStart'],
                                                  normalized['queryEnd'])
                e = hsp.expect
                if eCutoff is not None and e >= eCutoff:
                    continue

                if e == 0.0:
                    convertedE = None
                    plotInfo['zeroEValueFound'] = True
                else:
                    convertedE = -1.0 * log10(e)
                    if convertedE > plotInfo['maxE']:
                        plotInfo['maxE'] = convertedE
                    # Don't use elif for testing minE. Both can be true.
                    if convertedE < plotInfo['minE']:
                        plotInfo['minE'] = convertedE
                if normalized['queryStart'] < plotInfo['minX']:
                    plotInfo['minX'] = normalized['queryStart']
                if normalized['queryEnd'] > plotInfo['maxX']:
                    plotInfo['maxX'] = normalized['queryEnd']
                plotInfo['items'].append({
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

            if rankEValues:
                self._convertEValuesToRanks(plotInfo)

            # For each sequence we have hits on, set the expect values that
            # were zero to a randomly high value (higher than the max e value
            # we just calculated).
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

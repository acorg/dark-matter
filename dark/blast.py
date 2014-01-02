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
    @param fastaFilename: the C{str} file name containing the sequences that
        were given to BLAST as queries.
    """

    def __init__(self, blastFilename, fastaFilename):
        self._blastFilename = blastFilename
        self._fastaFilename = fastaFilename

    def records(self, limit=None):
        """
        Extract all BLAST records.

        @param limit: An C{int} limit on the number of records to return.

        @return: A generator that yields BioPython C{Bio.Blast.Record.Blast}
            instances.
        """
        if self._blastFilename.endswith('.xml'):
            fp = open(self._blastFilename)
            records = NCBIXML.parse(fp)
        elif self._blastFilename.endswith('.json'):
            records = readJSONRecords(self._blastFilename)
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
                        'length': alignment.length,  # The subject length.
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
            hitInfo['reads'] = sorted(hitInfo['reads'])
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
        self._records = records
        self._titles = {}
        self._fasta = None  # Computed (once) in summarizeHits.

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

        if title in self._titles:
            raise KeyError('Title %r already present' % (title,))

        self._titles[title] = hitInfo

    def titles(self):
        """
        Get the sequence titles these BLAST hits are for.

        @return: a generator yielding all hit titles.
        """
        return self._titles.iterkeys()

    def _sortBy(self, attr='count', reverse=False):
        """
        Return an HTML object with the hits sorted by the given attribute.

        @param attr: A C{str}, one of 'readCount', 'eMean', 'eMedian', or
            'length'. Raise C{ValueError} if C{attr} is not one of the above.
        @param reverse: If C{True}, reverse the sort order.
        """
        if attr not in ('readCount', 'eMean', 'eMedian', 'length'):
            raise ValueError("attr must be one of 'readCount', 'eMean', "
                             "'eMedian', or 'length'")

        out = []
        titles = sorted(self.titles(), reverse=reverse,
                        key=lambda title: self._titles[title][attr])
        for i, title in enumerate(titles, start=1):
            hitInfo = self._titles[title]
            link = NCBISequenceLink(title, title)
            out.append(
                '%3d: count=%4d, len=%7d, median(e)=%20s mean(e)=%20s: %s' %
                (i, hitInfo['readCount'], hitInfo['length'],
                 hitInfo['eMedian'], hitInfo['eMean'], link))
        return HTML('<pre><tt>' + '<br/>'.join(out) + '</tt></pre>')

    def summarizeByMeanEValueHTML(self):
        return self._sortBy('eMean')

    def summarizeByMedianEValueHTML(self):
        return self._sortBy('eMedian')

    def summarizeByCountHTML(self):
        return self._sortBy('readCount', reverse=True)

    def summarizeByLengthHTML(self):
        return self._sortBy('length', reverse=True)

    def filter(self, filterFunc):
        """
        Produce a new L{BlastHits} instance consisting of just the hits that
        C{filterFunc} returns a true value for.

        @param filterFunc: A C{function} that accepts two arguments, a string
        sequence title and a hitInfo dict (see docstring for C{addHit} above).
        The function should return a true value if the hit should appear in the
        result.

        NOTE: this function is not currently used (though there are tests for
              it). We should probably remove it.
        """
        result = BlastHits(self._records)
        for title, hitInfo in self._titles.iteritems():
            if filterFunc(title, hitInfo):
                result.addHit(title, hitInfo)
        return result

    def interesting(self, titleRegex=None, minSequenceLen=None,
                    maxSequenceLen=None, minMatchingReads=None,
                    maxMeanEValue=None, maxMedianEValue=None,
                    negativeTitleRegex=None, maxMinEValue=None,
                    truncateTitlesAfter=None):
        """
        Produce a new L{BlastHits} instance consisting of just the interesting
        hits, as given by our parameters.

        IMPORTANT: the hitInfo entries in the new instance are the same objects
        we have in self._titles. So if our caller changes the result we return,
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
        result = BlastHits(self._records)

        if truncateTitlesAfter:
            truncatedTitles = set()
        if titleRegex is not None:
            titleRegex = re.compile(titleRegex, re.I)
        if negativeTitleRegex is not None:
            negativeTitleRegex = re.compile(negativeTitleRegex, re.I)

        for title, hitInfo in self._titles.iteritems():
            if (minSequenceLen is not None and
                    hitInfo['length'] < minSequenceLen):
                continue
            if (maxSequenceLen is not None and
                    hitInfo['length'] > maxSequenceLen):
                continue
            if (minMatchingReads is not None and
                    hitInfo['count'] < minMatchingReads):
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
        titles = self._titles
        for readNum, record in enumerate(self._records.records()):
            for index, alignment in enumerate(record.alignments):
                title = record.descriptions[index].title
                if title in titles:
                    yield (title, readNum, alignment.hsps)

    def summarizeForPlotting(self, eCutoff=None, maxHspsPerHit=None,
                             minStart=None, maxStop=None, logLinearXAxis=False,
                             logBase=DEFAULT_LOG_LINEAR_X_AXIS_BASE,
                             randomizeZeroEValues=False):
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

        @return: A L{BlastHitsSummary} instance.
        """

        # Reset the plot info for each hit title because we may be called
        # multiple times, with different parameters.
        for title, hitInfo in self._titles.iteritems():
            hitInfo['plotInfo'] = None

        # Read in the read FASTA file if we haven't already.
        if self._fasta is None:
            self._fasta = list(SeqIO.parse(self._records._fastaFilename,
                                           'fasta'))

        zeroEValueUpperRandomIncrement = 150

        def plotInfoDict(sequenceLen):
            result = {
                'hspCount': 0,
                'items': [],
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

        for title, readNum, hsps in self._hitDetails():
            if self._titles[title]['plotInfo'] is None:
                self._titles[title]['plotInfo'] = plotInfoDict()
            plotInfo = self._titles[title]['plotInfo']
            queryLen = len(self._fasta[readNum])

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
                    # print 'hitLen: %s' % self._titles[title]['length']
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

                plotInfo['hspCount'] += 1
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

        for title in self._titles.iterkeys():
            plotInfo = self._titles[title]['plotInfo']
            if plotInfo is None:
                continue

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

                plotInfo.update({
                    'minX': minX,
                    'maxX': maxX,
                    'offsetAdjuster': adjuster,
                })


class BlastHitsSummary(object):
    """
    Hold summary information about a set of BLAST hits against a set of
    sequence titles.

    @param hits: A L{BlastHits} instance.
    """

    def __init__(self, hits):
        self._hits = hits

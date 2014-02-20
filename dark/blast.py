from math import log10
import numpy as np
from random import uniform
from Bio import SeqIO
import MySQLdb
from os import environ

from dark.conversion import JSONRecordsReader, XMLRecordsReader
from dark.filter import (BitScoreFilter, HitInfoFilter, ReadSetFilter,
                         TitleFilter)
from dark.hsp import printHSP, normalizeHSP
from dark.intervals import OffsetAdjuster, ReadIntervals

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

    def __init__(self, blastFilename, fastaFilename=None, blastDb=None,
                 limit=None):
        self.blastFilename = blastFilename
        self.fastaFilename = fastaFilename
        self.blastDb = blastDb
        self.limit = limit
        self._length = 0
        self.blastParams = None  # Set in records(), below.

    def records(self):
        """
        Extract all BLAST records (up to C{self.limit}, if not C{None}).

        @return: A generator that yields BioPython C{Bio.Blast.Record.Blast}
            instances.
        """
        if self.blastFilename.endswith('.xml'):
            reader = XMLRecordsReader(self.blastFilename)
        elif self.blastFilename.endswith('.json'):
            reader = JSONRecordsReader(self.blastFilename)
        else:
            raise ValueError('Unknown BLAST record file type.')

        # Read the records, observing any given limit. There's a little code
        # duplication here in the name of speed - we don't want to repeatedly
        # test if self.limit is None in a loop.
        count = 0
        if self.limit is None:
            for record in reader.records():
                count += 1
                yield record
        else:
            limit = self.limit
            for record in reader.records():
                if count == limit:
                    break
                count += 1
                yield record

        self.blastParams = reader.params
        self._length = count

    def __len__(self):
        """
        Return the number of BLAST records.

        NOTE: this will return zero if called before C{records()} has been
              called to read the records.

        @return: an C{int}, the number of records in this BLAST output. This is
            the number of reads that were present in the FASTA file given to
            BLAST because BLAST outputs a record for each input (query)
            sequence.
        """
        return self._length

    def filterHits(self, whitelist=None, blacklist=None, minSequenceLen=None,
                   maxSequenceLen=None, minMatchingReads=None,
                   maxMeanEValue=None, maxMedianEValue=None,
                   withEBetterThan=None, titleRegex=None,
                   negativeTitleRegex=None, truncateTitlesAfter=None,
                   minMeanBitScore=None, minMedianBitScore=None,
                   withBitScoreBetterThan=None, minNewReads=None):
        """
        Read the BLAST records and return a L{BlastHits} instance. Records are
        only returned if they match the various optional restrictions described
        below.

        @param whitelist: If not C{None}, a set of exact titles that are always
            acceptable (though the hit info for a whitelist title may rule it
            out for other reasons).
        @param blacklist: If not C{None}, a set of exact titles that are never
            acceptable.
        @param minSequenceLen: sequences of lesser length will be elided.
        @param maxSequenceLen: sequences of greater length will be elided.
        @param minMatchingReads: sequences that are matched by fewer reads
            will be elided.
        @param maxMeanEValue: sequences that are matched with a mean e-value
            that is greater will be elided.
        @param maxMedianEValue: sequences that are matched with a median
            e-value that is greater will be elided.
        @param withEBetterThan: if the best (minimum) e-value for a hit is not
            as good as (i.e., is higher than) this value, elide the hit. E.g.,
            suppose we are passed a value of 1e-20, then we should reject any
            hit whose best (i.e., lowest) e-value is worse (bigger) than 1e-20.
            So a hit with minimal e-value of 1e-10 would not be reported,
            whereas a hit with a minimal e-value of 1e-30 would be.
        @param titleRegex: a regex that sequence titles must match.
        @param negativeTitleRegex: a regex that sequence titles must not match.
        @param truncateTitlesAfter: specify a string that titles will be
            truncated beyond. If a truncated title has already been seen, that
            title will be elided.
        @param minMeanBitScore: sequences that are matched with a mean score
            that is less than this value will be elided.
        @param minMedianBitScore: sequences that are matched with a median
            score that is less than this value will be elided.
        @param withBitScoreBetterThan: If no score for a sequence is higher
            than this value, the hit will be elided.
        @param minNewReads: The C{float} fraction of its reads by which a new
            read set must differ from all previously seen read sets in order to
            be considered acceptably different.
        @return: A L{BlastHits} instance.
        """
        result = {}
        titleFilter = TitleFilter(
            whitelist=whitelist, blacklist=blacklist, positiveRegex=titleRegex,
            negativeRegex=negativeTitleRegex,
            truncateAfter=truncateTitlesAfter)

        # For each read (that BLAST found in the FASTA file)...
        for readNum, record in enumerate(self.records()):

            # For each sequence that the read matched against...
            for index, alignment in enumerate(record.alignments):

                # Test sequence title.
                title = record.descriptions[index].title
                titleFilterResult = titleFilter.accept(title)
                if titleFilterResult == TitleFilter.REJECT:
                    continue

                if title in result:
                    hitInfo = result[title]
                else:
                    # Test sequence length before adding it to result, to make
                    # sure the length tests are only done once per title.
                    sequenceLen = alignment.length
                    if ((minSequenceLen is not None and
                         sequenceLen < minSequenceLen) or
                        (maxSequenceLen is not None and
                         sequenceLen > maxSequenceLen)):
                        continue

                    hitInfo = result[title] = {
                        'eValues': [],
                        'length': sequenceLen,
                        'readCount': 0,
                        'readNums': set(),
                        'bitScores': [],
                        'titleFilterResult': titleFilterResult
                    }

                # Record just the best e-value with which this read hit this
                # sequence.
                hitInfo['eValues'].append(alignment.hsps[0].expect)
                hitInfo['bitScores'].append(alignment.hsps[0].bits)
                hitInfo['readCount'] += 1
                hitInfo['readNums'].add(readNum)

        # Note that we don't pass minSequenceLen or maxSequenceLen to the
        # hit info filter since we have already tested those.
        hitInfoFilter = HitInfoFilter(
            minMatchingReads=minMatchingReads, maxMeanEValue=maxMeanEValue,
            maxMedianEValue=maxMedianEValue, withEBetterThan=withEBetterThan)

        bitScoreFilter = BitScoreFilter(
            minMeanBitScore=minMeanBitScore,
            minMedianBitScore=minMedianBitScore,
            withBitScoreBetterThan=withBitScoreBetterThan)

        if minNewReads is None:
            readSetFilter = None
        else:
            readSetFilter = ReadSetFilter(minNewReads)

        # Compute summary stats on e-values for all titles. If the title
        # was whitelisted or if the statistical summary is acceptable, add
        # the hit info to our final result.

        blastHits = BlastHits(self, readSetFilter=readSetFilter)

        titles = result.keys()  # Don't change 'result' while we iterate it.
        for title in titles:
            hitInfo = result[title]
            eValues = hitInfo['eValues']
            hitInfo['eMean'] = np.mean(eValues)
            hitInfo['eMedian'] = np.median(eValues)
            hitInfo['eMin'] = min(eValues)
            bitScores = hitInfo['bitScores']
            hitInfo['bitScoreMean'] = np.mean(bitScores)
            hitInfo['bitScoreMedian'] = np.median(bitScores)
            hitInfo['bitScoreMax'] = np.max(bitScores)
            if (hitInfo['titleFilterResult'] == TitleFilter.WHITELIST_ACCEPT or
                    hitInfoFilter.accept(hitInfo) and
                    bitScoreFilter.accept(hitInfo) and
                    (minNewReads is None or
                     readSetFilter.accept(title, hitInfo))):
                # Remove the e-values and bit scores (now that we've summarized
                # them) and the title filter result (now that we've checked it.
                del hitInfo['eValues']
                del hitInfo['bitScores']
                del hitInfo['titleFilterResult']
                blastHits.addHit(title, hitInfo)

            # Reduce memory usage as quickly as we can.
            del result[title]

        return blastHits


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
        @param hitInfo: A C{dict} of information associated with hits (i.e.,
            matching reads) against the sequence with this title. The dict must
            contain the following keys:

                readCount: the number of reads that hit this sequence.
                length:    the length of the sequence.
                reads:     a sorted C{list} of read names (as found in e.g., a
                           FASTA file).
                eMean:     the mean of the e-values of all hits against this
                           sequence.
                eMedian:   the median of the e-values of all hits against this
                           sequence.
                eMin:      the minimum of the e-values of all hits against this
                           sequence.
        """

        if title in self.titles:
            raise KeyError('Title %r already present' % (title,))

        self.titles[title] = hitInfo

    def sortTitles(self, by):
        """
        Sort titles by a given attribute and then by title.

        @param by: A C{str}, one of 'eMean', 'eMedian', 'eMin', 'readCount',
            'title', 'length', 'bitScoreMax', 'bitScoreMedian', 'bitScoreMean'.
        @return: A sorted C{list} of titles.
        """

        def makeCmp(attr):
            """
            Create a sorting comparison function that sorts first in reverse
            on the passed numeric attribute and then in ascending order on
            title.
            """
            def compare(title1, title2):
                result = cmp(self.titles[title2][attr],
                             self.titles[title1][attr])
                if result == 0:
                    result = cmp(title1, title2)
                return result
            return compare

        if by == 'eMin':
            return sorted(
                self.titles.iterkeys(),
                key=lambda title: (self.titles[title]['eMin'], title))
        elif by == 'eMean':
            return sorted(
                self.titles.iterkeys(),
                key=lambda title: (self.titles[title]['eMean'], title))
        elif by == 'eMedian':
            return sorted(
                self.titles.iterkeys(),
                key=lambda title: (self.titles[title]['eMedian'], title))
        elif by == 'bitScoreMax':
            return sorted(
                self.titles.iterkeys(), cmp=makeCmp('bitScoreMax'))
        elif by == 'bitScoreMean':
            return sorted(
                self.titles.iterkeys(), cmp=makeCmp('bitScoreMean'))
        elif by == 'bitScoreMedian':
            return sorted(
                self.titles.iterkeys(), cmp=makeCmp('bitScoreMedian'))
        elif by == 'readCount':
            return sorted(self.titles.iterkeys(), cmp=makeCmp('readCount'))
        elif by == 'length':
            return sorted(self.titles.iterkeys(), cmp=makeCmp('length'))
        elif by == 'title':
            return sorted(self.titles.iterkeys())

        raise ValueError('sort attribute must be one of "eMean", '
                         '"eMedian", "eMin", "readCount", "title".')

    def sortTitlesOnPlotInfo(self, by):
        """
        Sort titles by a given attribute and then by title. To sort, use
        information in plotInfo.

        @param by: A C{str}, one of 'eMean', 'eMedian', 'eMin', 'readCount',
            'title', 'length', 'bitScoreMin', 'bitScoreMean', 'bitScoreMedian'.
        @return: A sorted C{list} of titles.
        """

        def makeCmp(attr):
            """
            Create a sorting comparison function that sorts first in reverse
            on the passed numeric attribute and then in ascending order on
            title.
            """
            def compare(title1, title2):
                result = cmp(self.titles[title2][attr],
                             self.titles[title1][attr])
                if result == 0:
                    result = cmp(title1, title2)
                return result
            return compare

        def makePlotInfoCmp(attr):
            """
            Create a sorting comparison function that sorts first in reverse
            on the passed numeric attribute and then in ascending order on
            title.
            """
            def compare(title1, title2):
                result = cmp(self.titles[title2]['plotInfo'][attr],
                             self.titles[title1]['plotInfo'][attr])
                if result == 0:
                    result = cmp(title1, title2)
                return result
            return compare

        def makePlotInfoKey(attr):
            """
            Create a function that returns a sorting key tuple consisting of
            the passed plotinfo attribute and then the sequence title.
            """
            def key(title):
                return (self.titles[title]['plotInfo'][attr], title)
            return key

        # Only return information about titles that have some plotinfo.
        titles = (title for title in self.titles
                  if self.titles[title]['plotInfo'] is not None)

        if by == 'eMin':
            return sorted(titles, key=makePlotInfoKey('originalEMin'))
        elif by == 'eMean':
            return sorted(titles, key=makePlotInfoKey('originalEMean'))
        elif by == 'eMedian':
            return sorted(titles, key=makePlotInfoKey('originalEMedian'))
        elif by == 'bitScoreMax':
            return sorted(titles, cmp=makePlotInfoCmp('bitScoreMax'))
        elif by == 'bitScoreMean':
            return sorted(titles, cmp=makePlotInfoCmp('bitScoreMean'))
        elif by == 'bitScoreMedian':
            return sorted(titles, cmp=makePlotInfoCmp('bitScoreMedian'))
        elif by == 'readCount':
            return sorted(titles, cmp=makeCmp('readCount'))
        elif by == 'length':
            return sorted(titles, cmp=makeCmp('length'))
        elif by == 'title':
            return sorted(titles)

        raise ValueError('sort attribute must be one of "eMean", '
                         '"eMedian", "eMin", "bitScoreMax", "bitScoreMean", '
                         '"bitScoreMedian", "readCount", "length", "title".')

    def filterHits(self, whitelist=None, blacklist=None, minSequenceLen=None,
                   maxSequenceLen=None, minMatchingReads=None,
                   maxMeanEValue=None, maxMedianEValue=None,
                   withEBetterThan=None, titleRegex=None,
                   negativeTitleRegex=None, truncateTitlesAfter=None,
                   minMeanBitScore=None, minMedianBitScore=None,
                   withBitScoreBetterThan=None, minNewReads=None):
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

        @param whitelist: If not C{None}, a set of exact titles that are always
            acceptable (though the hit info for a whitelist title may rule it
            out for other reasons).
        @param blacklist: If not C{None}, a set of exact titles that are never
            acceptable.
        @param minSequenceLen: sequences of lesser length will be elided.
        @param maxSequenceLen: sequences of greater length will be elided.
        @param minMatchingReads: sequences that are matched by fewer reads
            will be elided.
        @param maxMeanEValue: sequences that are matched with a mean e-value
            that is greater will be elided.
        @param maxMedianEValue: sequences that are matched with a median
            e-value that is greater will be elided.
        @param withEBetterThan: if the best (minimum) e-value for a hit is not
            as good as (i.e., is higher than) this value, elide the hit. E.g.,
            suppose we are passed a value of 1e-20, then we should reject any
            hit whose best (i.e., lowest) e-value is worse (bigger) than 1e-20.
            So a hit with minimal e-value of 1e-10 would not be reported,
            whereas a hit with a minimal e-value of 1e-30 would be.
        @param titleRegex: a regex that sequence titles must match.
        @param negativeTitleRegex: a regex that sequence titles must not match.
        @param truncateTitlesAfter: specify a string that titles will be
            truncated beyond. If a truncated title has already been seen, that
            title will be elided.
        @param minMeanBitScore: sequences that are matched with a mean bit
            score that is less than this value will be elided.
        @param minMedianBitScore: sequences that are matched with a median
            bit score that is less than this value will be elided.
        @param withBitScoreBetterThan: If no bit score for a sequence is higher
            than this value, the hit will be elided.
        @param minNewReads: The C{float} fraction of its reads by which a new
            read set must differ from all previously seen read sets in order to
            be considered acceptably different.
        @return: A new L{BlastHits} instance, with hits filtered as above.
        """
        titleFilter = TitleFilter(
            whitelist=whitelist, blacklist=blacklist, positiveRegex=titleRegex,
            negativeRegex=negativeTitleRegex,
            truncateAfter=truncateTitlesAfter)

        hitInfoFilter = HitInfoFilter(
            minSequenceLen=minSequenceLen, maxSequenceLen=maxSequenceLen,
            minMatchingReads=minMatchingReads, maxMeanEValue=maxMeanEValue,
            maxMedianEValue=maxMedianEValue, withEBetterThan=withEBetterThan)

        bitScoreFilter = BitScoreFilter(
            minMeanBitScore=minMeanBitScore,
            minMedianBitScore=minMedianBitScore,
            withBitScoreBetterThan=withBitScoreBetterThan)

        # Use a ReadSetFilter only if we're checking that read sets are
        # sufficiently new.
        if minNewReads is None:
            readSetFilter = None
        else:
            readSetFilter = ReadSetFilter(minNewReads)

        blastHits = BlastHits(self.records, readSetFilter=readSetFilter)
        for title, hitInfo in self.titles.iteritems():
            titleFilterResult = titleFilter.accept(title)
            if (titleFilterResult == TitleFilter.WHITELIST_ACCEPT or
                    titleFilterResult == TitleFilter.DEFAULT_ACCEPT and
                    hitInfoFilter.accept(hitInfo) and
                    bitScoreFilter.accept(hitInfo) and
                    (minNewReads is None or
                     readSetFilter.accept(title, hitInfo))):
                blastHits.addHit(title, hitInfo)
        return blastHits

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
                if title in titles and readNum in titles[title]['readNums']:
                    yield (title, readNum, alignment.hsps)

    def _convertEValuesToRanks(self, plotInfo):
        """
        Change e-values for the reads that hit each sequence to be their ranks.

        @param plotInfo: A C{dict} of plot information, as returned by
            C{plotInfoDict} in C{computePlotInfo}.
        """
        items = plotInfo['items']
        items.sort(key=lambda item: item['convertedE'])
        for rank, item in enumerate(items, start=1):
            item['convertedE'] = rank
        plotInfo['maxE'] = plotInfo['maxEIncludingRandoms'] = len(items)
        plotInfo['minE'] = 1
        plotInfo['zeroEValueFound'] = False

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
                        randomizeZeroEValues=False, rankValues=False):
        """
        Read detailed HSP information about the hit titles in C{self.titles}
        and compute summary statistics on it. The various parameters allow
        us to restrict and transform the data that is read.

        @param eCutoff: A float e-value. Hits with e-value greater than or
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
            'eCutoff': eCutoff,
            'logBase': logBase,
            'logLinearXAxis': logLinearXAxis,
            'maxHspsPerHit': maxHspsPerHit,
            'maxStop': maxStop,
            'minStart': minStart,
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
                # these two below should go away
                'eMean': None,
                'eMedian': None,
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

            if maxHspsPerHit is not None and len(hsps) > maxHspsPerHit:
                hsps = hsps[:maxHspsPerHit]

            for hspCount, hsp in enumerate(hsps, start=1):
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

                # We are committed to adding this item to plotInfo['items'].
                # Don't add any 'continue' statements below this point.

                # retain original evalue below cutoff for eMean and eMedian
                # calculation.
                plotInfo['originalEValues'].append(e)

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

    def getTaxIDFromMySql(self):
        """
        For each title in C{self.titles}, read the corresponding taxId from
        the gi_taxid_nucl table in a MySQL database called ncbi_taxonomy.
        """
        # connect to database
        # parameters should be changed accordingly
        db = MySQLdb.connect(host='localhost', user=environ.get(
                             'DBI_USER', environ['USER']),
                             passwd=environ['DBI_PASSWORD'],
                             db='ncbi_taxonomy')
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
            if result == 'No taxID found':
                print title
            self.titles[title]['taxID'] = result

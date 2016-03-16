from __future__ import print_function

import bz2
from json import dumps, loads
from operator import itemgetter

from Bio.Blast import NCBIXML
from Bio.File import as_handle

from dark.hsp import HSP, LSP
from dark.score import HigherIsBetterScore
from dark.alignments import Alignment, ReadAlignments
from dark.blast.hsp import normalizeHSP


class XMLRecordsReader(object):
    """
    Provide a method that yields parsed XML records from a file. Store and
    make accessible the global BLAST parameters.

    @ivar params: A C{dict} of global BLAST parameters.
    @param filename: A C{str} filename or an open file pointer, containing XML
        BLAST records.
    """

    def __init__(self, filename):
        self._filename = filename
        self.params = None  # Set below, in records.

    def _convertBlastRecordToDict(self, record):
        """
        Pull (only) the fields we use out of the record and return them as a
        dict.  Although we take the title from each alignment description, we
        save space in the JSON output by storing it in the alignment dict (not
        in a separated 'description' dict). When we undo this conversion (in
        JSONRecordsReader._convertDictToBlastRecord) we'll pull the title out
        of the alignment dict and put it into the right place in the BLAST
        record.

        @param record: An instance of C{Bio.Blast.Record.Blast}. The attributes
            on this don't seem to be documented. You'll need to look at the
            BioPython source to see everything it contains.
        @return: A C{dict} with 'alignments' and 'query' keys.
        """
        alignments = []
        for alignment in record.alignments:
            hsps = []
            for hsp in alignment.hsps:
                hsps.append({
                    'bits': hsp.bits,
                    'expect': hsp.expect,
                    'frame': hsp.frame,
                    'query': hsp.query,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'sbjct': hsp.sbjct,
                    'sbjct_start': hsp.sbjct_start,
                    'sbjct_end': hsp.sbjct_end,
                })

            alignments.append({
                'hsps': hsps,
                'length': alignment.length,
                'title': alignment.title,
            })

        return {
            'alignments': alignments,
            'query': record.query,
        }

    def _convertBlastParamsToDict(self, record):
        """
        Pull the global BLAST parameters out of a BLAST record and return
        them as a C{dict}.

        Some of these attributes are useless (not filled in), but we record
        them all just in case we one day need them or they start to be used or
        they disappear etc. Any of those changes might alert us that something
        has changed in BLAST XML output or in BioPython.

        @param record: An instance of C{Bio.Blast.Record.Blast}. The attributes
            on this don't seem to be documented. You'll need to look at the
            BioPython source to see everything it contains.
        @return: A C{dict}, as described above.
        """
        result = {}
        for attr in (
                # From Bio.Blast.Record.Header
                'application',
                'version',
                'date',
                'reference',
                'query',
                'query_letters',
                'database',
                'database_sequences',
                'database_letters',
                # From Bio.Blast.Record.DatabaseReport
                'database_name',
                'posted_date',
                'num_letters_in_database',
                'num_sequences_in_database',
                'ka_params',
                'gapped',
                'ka_params_gap',
                # From Bio.Blast.Record.Parameters
                'matrix',
                'gap_penalties',
                'sc_match',
                'sc_mismatch',
                'num_hits',
                'num_sequences',
                'num_good_extends',
                'num_seqs_better_e',
                'hsps_no_gap',
                'hsps_prelim_gapped',
                'hsps_prelim_gapped_attemped',
                'hsps_gapped',
                'query_id',
                'query_length',
                'database_length',
                'effective_hsp_length',
                'effective_query_length',
                'effective_database_length',
                'effective_search_space',
                'effective_search_space_used',
                'frameshift',
                'threshold',
                'window_size',
                'dropoff_1st_pass',
                'gap_x_dropoff',
                'gap_x_dropoff_final',
                'gap_trigger',
                'blast_cutoff'):
            result[attr] = getattr(record, attr)
        return result

    def records(self):
        """
        Yield BLAST records, as read by the BioPython NCBIXML.parse
        method. Set self.params from data in the first record.
        """
        first = True
        with as_handle(self._filename) as fp:
            for record in NCBIXML.parse(fp):
                if first:
                    self.params = self._convertBlastParamsToDict(record)
                    first = False
                yield record

    def saveAsJSON(self, fp):
        """
        Write the records out as JSON. The first JSON object saved contains
        the BLAST parameters.

        @param fp: A C{str} file pointer to write to.
        """
        first = True
        for record in self.records():
            if first:
                print(dumps(self.params, separators=(',', ':')), file=fp)
                first = False
            print(dumps(self._convertBlastRecordToDict(record),
                        separators=(',', ':')), file=fp)


class JSONRecordsReader(object):
    """
    Provide a method that yields JSON records from a file. Store, check, and
    make accessible the global BLAST parameters.

    @param filename: A C{str} filename containing JSON BLAST records.
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}, for comparing bit scores. If you
        are using e-values, pass LowerIsBetterScore instead.
    """

    # Note that self._fp is opened in self.__init__, accessed in
    # self._params and in self.records, and closed in self.close.

    def __init__(self, filename, scoreClass=HigherIsBetterScore):
        self._filename = filename
        self._scoreClass = scoreClass
        if scoreClass is HigherIsBetterScore:
            self._hspClass = HSP
        else:
            self._hspClass = LSP

        self._open(filename)
        self.application = self.params['application'].lower()

    def _open(self, filename):
        """
        Open the input file. Set self._fp to point to it. Read the first
        line of parameters.

        @param filename: A C{str} filename containing JSON BLAST records.
        @raise ValueError: if the first line of the file isn't valid JSON,
            if the input file is empty, or if the JSON does not contain an
            'application' key.
        """
        if filename.endswith('.bz2'):
            self._fp = bz2.BZ2File(filename)
        else:
            self._fp = open(filename)

        line = self._fp.readline()
        if not line:
            raise ValueError('JSON file %r was empty.' % self._filename)

        try:
            self.params = loads(line[:-1])
        except ValueError as e:
            raise ValueError(
                'Could not convert first line of %r to JSON (%s). '
                'Line is %r.' % (self._filename, e, line[:-1]))
        else:
            if 'application' not in self.params:
                raise ValueError(
                    '%r appears to be an old JSON file with no BLAST global '
                    'parameters. Please re-run convert-blast-xml-to-json.py '
                    'to convert it to the newest format.' % self._filename)

    def _dictToAlignments(self, blastDict, read):
        """
        Take a dict (made by XMLRecordsReader._convertBlastRecordToDict)
        and convert it to a list of alignments.

        @param blastDict: A C{dict}, from convertBlastRecordToDict.
        @param read: A C{Read} instance, containing the read that BLAST used
            to create this record.
        @raise ValueError: If the query id in the BLAST dictionary does not
            match the id of the read.
        @return: A C{list} of L{dark.alignment.Alignment} instances.
        """
        if (blastDict['query'] != read.id and
                blastDict['query'].split()[0] != read.id):
            raise ValueError(
                'The reads you have provided do not match the BLAST output: '
                'BLAST record query id (%s) does not match the id of the '
                'supposedly corresponding read (%s).' %
                (blastDict['query'], read.id))

        alignments = []
        getScore = itemgetter('bits' if self._hspClass is HSP else 'expect')

        for blastAlignment in blastDict['alignments']:
            alignment = Alignment(blastAlignment['length'],
                                  blastAlignment['title'])
            alignments.append(alignment)
            for blastHsp in blastAlignment['hsps']:
                score = getScore(blastHsp)
                normalized = normalizeHSP(blastHsp, len(read),
                                          self.application)
                hsp = self._hspClass(
                    score,
                    readStart=normalized['readStart'],
                    readEnd=normalized['readEnd'],
                    readStartInSubject=normalized['readStartInSubject'],
                    readEndInSubject=normalized['readEndInSubject'],
                    readFrame=blastHsp['frame'][0],
                    subjectStart=normalized['subjectStart'],
                    subjectEnd=normalized['subjectEnd'],
                    subjectFrame=blastHsp['frame'][1],
                    readMatchedSequence=blastHsp['query'],
                    subjectMatchedSequence=blastHsp['sbjct'])

                alignment.addHsp(hsp)

        return alignments

    def readAlignments(self, reads):
        """
        Read lines of JSON from self._filename, convert them to read alignments
        and yield them.

        @param reads: An iterable of L{Read} instances, corresponding to the
            reads that were given to BLAST.
        @raise ValueError: If any of the lines in the file cannot be converted
            to JSON.
        @return: A generator that yields C{dark.alignments.ReadAlignments}
            instances.
        """
        if self._fp is None:
            self._open(self._filename)

        reads = iter(reads)

        try:
            for lineNumber, line in enumerate(self._fp, start=2):
                try:
                    record = loads(line[:-1])
                except ValueError as e:
                    raise ValueError(
                        'Could not convert line %d of %r to JSON (%s). '
                        'Line is %r.' %
                        (lineNumber, self._filename, e, line[:-1]))
                else:
                    try:
                        read = next(reads)
                    except StopIteration:
                        raise ValueError(
                            'Read generator failed to yield read number %d '
                            'during parsing of BLAST file %r.' %
                            (lineNumber - 1, self._filename))
                    else:
                        alignments = self._dictToAlignments(record, read)
                        yield ReadAlignments(read, alignments)
        finally:
            self._fp.close()
            self._fp = None

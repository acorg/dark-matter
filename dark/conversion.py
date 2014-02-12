from Bio.Blast import NCBIXML
from Bio.Blast.Record import Alignment, Blast, Description, HSP
from json import dumps, loads


class XMLRecordsReader(object):
    """
    Provide a method that yields parsed XML records from a file. Store and
    make accessible the global BLAST parameters.

    @ivar params: A C{dict} of global BLAST parameters.
    @param filename: A C{str} filename containing XML BLAST records.
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
        for index, alignment in enumerate(record.alignments):
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
                'title': record.descriptions[index].title,
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
                ## From Bio.Blast.Record.Header
                'application',
                'version',
                'date',
                'reference',
                'query',
                'query_letters',
                'database',
                'database_sequences',
                'database_letters',
                ## From Bio.Blast.Record.DatabaseReport
                'database_name',
                'posted_date',
                'num_letters_in_database',
                'num_sequences_in_database',
                'ka_params',
                'gapped',
                'ka_params_gap',
                ## From Bio.Blast.Record.Parameters
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
        with open(self._filename) as fp:
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
                print >>fp, dumps(self.params, separators=(',', ':'))
                first = False
            print >>fp, dumps(self._convertBlastRecordToDict(record),
                              separators=(',', ':'))


class JSONRecordsReader(object):
    """
    Provide a method that yields JSON records from a file. Store, check, and
    make accessible the global BLAST parameters.

    @ivar params: A C{dict} of global BLAST parameters.
    @param filename: A C{str} filename containing JSON BLAST records.
    """

    def __init__(self, filename):
        self._filename = filename
        self.params = None  # Set below, in records.

    def _convertDictToBlastRecord(self, d):
        """
        Take a dictionary (as produced by
        XMLRecordsReader._convertBlastRecordToDict) and convert it to a Bio
        Blast record.

        @param d: A C{dict}, from convertBlastRecordToDict.
        @return: A C{Bio.Blast.Record.Blast} instance.
        """
        record = Blast()
        record.query = d['query']

        for alignment in d['alignments']:
            alignmentInstance = Alignment()
            for hsp in alignment['hsps']:
                hspInstance = HSP()
                for attr in ['bits', 'expect', 'frame', 'query', 'query_start',
                             'query_end', 'sbjct', 'sbjct_start', 'sbjct_end']:
                    setattr(hspInstance, attr, hsp[attr])
                alignmentInstance.hsps.append(hspInstance)
            title = alignment['title']
            alignmentInstance.hit_id = title.split(' ', 1)[0]
            alignmentInstance.length = alignment['length']
            record.alignments.append(alignmentInstance)

            descriptionInstance = Description()
            descriptionInstance.title = title
            record.descriptions.append(descriptionInstance)

        return record

    def _processFirstLine(self, line):
        """
        Look for BLAST parameters in the first line of an input file.

        @param line: A C{str} line of input, with a trailing '\n'. The
            line is expected to encode a JSON object.
        @raise ValueError: If the input file is empty or does not contain a
            valid JSON BLAST params section.
        @return: a C{dict} corresponding to the parsed JSON input line.
        """
        if not line:
            raise ValueError('Input file %r was empty' % self._filename)

        try:
            params = loads(line[:-1])
        except ValueError as e:
            raise ValueError(
                'Could not convert first line of %r to JSON (%s). Line is '
                '%r.' % (self._filename, e, line[:-1]))
        else:
            if 'application' not in params:
                raise ValueError(
                    '%r appears to be an old JSON file with no BLAST global '
                    'parameters. Please re-run convert-blast-xml-to-json.py '
                    'to convert it to the newest format.' % self._filename)
            return params

    def records(self):
        """
        Read lines of JSON from self._filename, convert them to Bio Blast class
        instances and yield them.

        @raise ValueError: If the input file is empty or does not contain a
            global BLAST params section.
        """
        with open(self._filename) as fp:
            for lineNumber, line in enumerate(fp, start=1):
                if lineNumber == 1:
                    self.params = self._processFirstLine(line)
                else:
                    try:
                        record = loads(line[:-1])
                    except ValueError as e:
                        raise ValueError(
                            'Could not convert line %d of %r to JSON (%s). '
                            'Line is %r.' %
                            (lineNumber, self._filename, e, line[:-1]))
                    else:
                        if 'application' in record:
                            # This is another params section. Check it
                            # against the already read params (from the
                            # first line of the file).  The 'application'
                            # key is just one of the many attributes we
                            # save into a JSON dict in
                            # L{XMLRecordsReader.convertBlastParamsToDict}
                            # above.
                            #
                            # Note that although the params contains a
                            # 'date', its value is empty (as far as I've
                            # seen). This could become an issue one day if
                            # it becomes non-empty and differs between JSON
                            # files that we cat together. In that case we
                            # may need to be more specific in our params
                            # compatible checking.
                            #
                            # TODO: We can be more helpful here by
                            #       reporting the keys that differ between
                            #       the params sections.
                            assert self.params == record, (
                                'BLAST parameters found on line %d of %r do '
                                'not match those found on its first line.' %
                                (lineNumber, self._filename))
                        else:
                            # A regular BLAST record (as a dict).
                            yield self._convertDictToBlastRecord(record)

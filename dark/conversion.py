from Bio.Blast import NCBIXML
from Bio.Blast.Record import Alignment, Blast, Description, HSP
from json import dumps, loads


def convertBlastXMLToJSON(blastFilename, jsonFilename):
    """
    Read BLAST XML records from blastFilename and write them
    out as serialized JSON to jsonFilename.
    """
    def _write(infp, outfp):
        first = True
        for record in NCBIXML.parse(infp):
            if first:
                print >>outfp, dumps(convertBlastParamsToDict(record))
                first = False
            print >>outfp, dumps(convertBlastRecordToDict(record))

    with open(blastFilename) as infp:
        if isinstance(jsonFilename, str):
            with open(jsonFilename, 'w') as outfp:
                _write(infp, outfp)
        else:
            # jsonFilename is already a file pointer (e.g., sys.stdout).
            _write(infp, jsonFilename)


def convertBlastParamsToDict(record):
    """
    Pull the global BLAST parameters out of a BLAST record and return
    them as a C{dict}.

    Some of these attributes are useless (not filled in), but we record them
    all just in case we one day need them or they start to be used or they
    disappear etc. Any of those changes might alert us that something has
    changed in BLAST XML output or in BioPython.

    @param record: An instance of C{Bio.Blast.Record.Blast}. The attributes
        on this don't seem to be documented. You'll need to look at the
        BioPython source to see everything it contains.
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


def convertBlastRecordToDict(record):
    """
    Pull (only) the fields we use out of the record and return them as a
    dict.  Although we take the title from each alignment description, we
    save space in the JSON output by storing it in the alignment dict (not
    in a separated 'description' dict). When we undo this conversion (in
    convertDictToBlastRecord) we'll pull the title out of the alignment
    dict and put it into the right place in the BLAST record.
    """
    alignments = []
    for index, alignment in enumerate(record.alignments):
        hsps = []
        for hsp in alignment.hsps:
            hsps.append({
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


def convertDictToBlastRecord(d):
    """
    Take a dictionary (as produced by convertBlastRecordToDict) and
    convert it to a Bio Blast record.
    """
    record = Blast()
    record.query = d['query']

    for alignment in d['alignments']:
        alignmentInstance = Alignment()
        for hsp in alignment['hsps']:
            hspInstance = HSP()
            for attr in ['expect', 'frame', 'query', 'query_start',
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


def readJSONRecords(filename):
    """
    Read lines of JSON from filename, convert them to Bio Blast class
    instances and yield them.
    """
    with open(filename) as fp:
        for lineNumber, line in enumerate(fp.readlines(), start=1):
            try:
                record = loads(line[:-1])
            except ValueError as e:
                raise ValueError(
                    'Could not convert line %d of %r to JSON (%s). Line is '
                    '%r.' % (lineNumber, filename, e, line[:-1]))
            else:
                yield convertDictToBlastRecord(record)

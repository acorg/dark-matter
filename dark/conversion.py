from Bio.Blast import NCBIXML
from Bio.Blast.Record import Alignment, Blast, Description, HSP
from json import dumps, loads


def convertBlastXMLToJSON(blastFilename, jsonFilename):
    """
    Read BLAST XML records from blastFilename and write them
    out as serialized JSON to jsonFilename.
    """
    def _write(infp, outfp):
        for record in NCBIXML.parse(infp):
            print >>outfp, dumps(convertBlastRecordToDict(record))

    with open(blastFilename) as infp:
        if isinstance(jsonFilename, str):
            with open(jsonFilename, 'w') as outfp:
                _write(infp, outfp)
        else:
            # jsonFilename is already a file pointer (e.g., sys.stdout).
            _write(infp, jsonFilename)


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
            except ValueError:
                raise ValueError(
                    'Could not convert line %d of %r to JSON. Line is %r.' %
                    (lineNumber, filename, line[:-1]))
            else:
                yield convertDictToBlastRecord(record)

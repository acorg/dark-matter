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
    dict.
    """
    descriptions = []
    for description in record.descriptions:
        descriptions.append({
            'e': description.e,
            'title': description.title,
        })

    alignments = []
    for alignment in record.alignments:
        hsps = []
        for hsp in alignment.hsps:
            hsps.append({
                'expect': hsp.expect,
                'frame': hsp.frame,
                'query': hsp.query,
                'query_start': hsp.query_start,
                'query_end': hsp.query_end,
                'match': hsp.match,
                'sbjct': hsp.sbjct,
                'sbjct_start': hsp.sbjct_start,
                'sbjct_end': hsp.sbjct_end,
            })

        alignments.append({
            'hit_id': alignment.hit_id,
            'hsps': hsps,
            'length': alignment.length,
        })

    return {
        'alignments': alignments,
        'descriptions': descriptions,
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
                         'query_end', 'match', 'sbjct', 'sbjct_start',
                         'sbjct_end']:
                setattr(hspInstance, attr, hsp[attr])
            alignmentInstance.hsps.append(hspInstance)
        alignmentInstance.hit_id = alignment['hit_id']
        alignmentInstance.length = alignment['length']
        record.alignments.append(alignmentInstance)

    for description in d['descriptions']:
        descriptionInstance = Description()
        descriptionInstance.e = description['e']
        descriptionInstance.title = description['title']
        record.descriptions.append(descriptionInstance)

    return record


def readJSONRecords(filename):
    """
    Read lines of JSON from filename, convert them to Bio Blast class
    instances and yield them.
    """
    for line in open(filename).readlines():
        record = loads(line[:-1])
        yield convertDictToBlastRecord(record)

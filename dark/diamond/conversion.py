from __future__ import print_function

import six
import bz2
from json import dumps, loads
from operator import itemgetter

from Bio.File import as_handle

from dark.hsp import HSP, LSP
from dark.score import HigherIsBetterScore
from dark.alignments import Alignment, ReadAlignments
from dark.diamond.hsp import normalizeHSP


class DiamondTabularFormatReader(object):
    """
    Provide a method that yields parsed tabular records from a file. Store and
    make accessible the global DIAMOND parameters.

    Make sure you run DIAMOND with the right output format. You must use:
        --outfmt "6 qtitle stitle bitscore evalue qframe qseq qstart qend sseq
                    sstart send slen"

    @param filename: A C{str} filename or an open file pointer, containing
        DIAMOND tabular records.
    """

    def __init__(self, filename):
        self._filename = filename
        self.application = 'DIAMOND'
        self.params = {
            'reference': ('Buchfink et al., Fast and Sensitive Protein '
                          'Alignment using DIAMOND, Nature Methods, 12, 59-60 '
                          '(2015)'),
            'task': 'blastx',  # TODO: Add support for blastp, if needed.
            'version': 'v0.8.23',
        }

    def records(self):
        """
        Parse the DIAMOND output and yield records.

        @return: A generator that produces C{dict}s containing 'alignments' and
            'query' C{str} keys.
        """
        with as_handle(self._filename) as fp:
            previousQseqid = None
            subjectsSeen = {}
            record = {}
            for line in fp:
                (qseqid, sseqid, bitscore, evalue, qframe, qseq, qstart, qend,
                 sseq, sstart, send, slen) = line.split('\t')
                hsp = {
                    'bits': float(bitscore),
                    'expect': float(evalue),
                    'frame': int(qframe),
                    'query': qseq,
                    'query_start': int(qstart),
                    'query_end': int(qend),
                    'sbjct': sseq,
                    'sbjct_start': int(sstart),
                    'sbjct_end': int(send),
                }
                if previousQseqid == qseqid:
                    # We have already started accumulating alignments for this
                    # query.
                    if sseqid not in subjectsSeen:
                        # We have not seen this subject before, so this is a
                        # new alignment.
                        subjectsSeen.add(sseqid)
                        alignment = {
                            'hsps': [hsp],
                            'length': int(slen),
                            'title': sseqid,
                        }
                        record['alignments'].append(alignment)
                    else:
                        # We have already seen this subject, so this is another
                        # HSP in an already existing alignment.
                        for alignment in record['alignments']:
                            if alignment['title'] == sseqid:
                                alignment['hsps'].append(hsp)
                                break
                else:
                    # All alignments for the previous query id (if any)
                    # have been seen.
                    if previousQseqid is not None:
                        yield record

                    # Start building up the new record.
                    record = {}
                    subjectsSeen = {sseqid}
                    alignment = {
                        'hsps': [hsp],
                        'length': int(slen),
                        'title': sseqid,
                    }
                    record['alignments'] = [alignment]
                    record['query'] = qseqid

                    previousQseqid = qseqid

            # Yield the last record, if any.
            if record:
                yield record

    def saveAsJSON(self, fp, writeBytes=False):
        """
        Write the records out as JSON. The first JSON object saved contains
        information about the DIAMOND algorithm.

        @param fp: A C{str} file pointer to write to.
        @param writeBytes: If C{True}, the JSON will be written out as bytes
            (not strings). This is required when fp is a bz2file.BZ2File.
        """
        if writeBytes:
            fp.write(dumps(self.params, sort_keys=True).encode('UTF-8'))
            fp.write(b'\n')
            for record in self.records():
                fp.write(dumps(record, sort_keys=True).encode('UTF-8'))
                fp.write(b'\n')
        else:
            fp.write(six.u(dumps(self.params, sort_keys=True)))
            fp.write(six.u('\n'))
            for record in self.records():
                fp.write(six.u(dumps(record, sort_keys=True)))
                fp.write(six.u('\n'))


class JSONRecordsReader(object):
    """
    Provide a method that yields JSON records from a file. Store, check, and
    make accessible the DIAMOND parameters.

    @param filename: A C{str} filename containing JSON DIAMOND records.
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}, for comparing bit scores. If you
        are using e-values, pass LowerIsBetterScore instead.
    """
    def __init__(self, filename, scoreClass=HigherIsBetterScore):
        self._filename = filename
        self._scoreClass = scoreClass
        if scoreClass is HigherIsBetterScore:
            self._hspClass = HSP
        else:
            self._hspClass = LSP

        self._open(filename)
        self.diamondTask = self.params['task']

    def _open(self, filename):
        """
        Open the input file. Set self._fp to point to it. Read the first
        line of parameters.

        @param filename: A C{str} filename containing JSON DIAMOND records.
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
                    '%r appears to be an old JSON file with no DIAMOND global '
                    'parameters. Please re-run convert-diamond-to-json.py '
                    'to convert it to the newest format.' % self._filename)

    def _dictToAlignments(self, diamondDict, read):
        """
        Take a dict (made by DiamondTabularFormatReader.records)
        and convert it to a list of alignments.

        @param diamondDict: A C{dict}, from records().
        @param read: A C{Read} instance, containing the read that DIAMOND used
            to create this record.
        @raise ValueError: If the query id in the DIAMOND dictionary does not
            match the id of the read.
        @return: A C{list} of L{dark.alignment.Alignment} instances.
        """
        if (diamondDict['query'] != read.id and
                diamondDict['query'].split()[0] != read.id):
            raise ValueError(
                'The reads you have provided do not match the DIAMOND output: '
                'DIAMOND record query id (%s) does not match the id of the '
                'supposedly corresponding read (%s).' %
                (diamondDict['query'], read.id))

        alignments = []
        getScore = itemgetter('bits' if self._hspClass is HSP else 'expect')

        for diamondAlignment in diamondDict['alignments']:
            alignment = Alignment(diamondAlignment['length'],
                                  diamondAlignment['title'])
            alignments.append(alignment)
            for diamondHsp in diamondAlignment['hsps']:
                score = getScore(diamondHsp)
                normalized = normalizeHSP(diamondHsp, len(read),
                                          self.diamondTask)
                hsp = self._hspClass(
                    score,
                    readStart=normalized['readStart'],
                    readEnd=normalized['readEnd'],
                    readStartInSubject=normalized['readStartInSubject'],
                    readEndInSubject=normalized['readEndInSubject'],
                    readFrame=diamondHsp['frame'],
                    subjectStart=normalized['subjectStart'],
                    subjectEnd=normalized['subjectEnd'],
                    readMatchedSequence=diamondHsp['query'],
                    subjectMatchedSequence=diamondHsp['sbjct'])

                alignment.addHsp(hsp)

        return alignments

    def readAlignments(self, reads):
        """
        Read lines of JSON from self._filename, convert them to read alignments
        and yield them.

        @param reads: An iterable of L{Read} instances, corresponding to the
            reads that were given to DIAMOND.
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
                    while True:
                        # Iterate through the input reads until we find the
                        # one that matches this DIAMOND record.
                        try:
                            read = next(reads)
                        except StopIteration:
                            raise ValueError(
                                'Read generator failed to yield a read '
                                'with id %r as found in record number %d '
                                'during parsing of DIAMOND output file %r.' %
                                (record['query'], lineNumber - 1,
                                 self._filename))
                        else:
                            if read.id == record['query']:
                                alignments = self._dictToAlignments(record,
                                                                    read)
                                yield ReadAlignments(read, alignments)
                                break
                            else:
                                # This is an input read that received no
                                # matches from DIAMOND. So it does not
                                # appear in the DIAMOND output. Emit an
                                # empty ReadAlignments for it.
                                yield ReadAlignments(read, [])

        finally:
            self._fp.close()
            self._fp = None

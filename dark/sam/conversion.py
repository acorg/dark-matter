from json import dumps
from itertools import groupby
import explain_sam_flags
import subprocess as sp
from collections import defaultdict
from dark.alignments import (Alignment, ReadAlignments)
from dark.hsp import HSP


class SAMtoJSON(object):  # Kept this in just in case, not debugged.
    """
    Provide a method that yields parsed SAM records from a file. Store as
    a JSON file.

    @param filename: A C{str} filename containing SAM records.
    @param outFile: A C{str} filename of a JSON file.
    """

    def __init__(self, filename, outFile):
        self.filename = filename
        self.outFile = outFile
        self.params = None  # Set below, in records.

    def _convertCigarMD(self, category, cigar, MD=None, match=6,
                        mismatch=-3, insert=-1, delete=-1):
        """
        Calculates sequence identity or 'bitscore' from CIGAR
        and MD string.
        Bitscore scoring: =:6, X:-3, I:-1, D:-1, can be changed.

        @param category: a C{str} of either 'seqID' or 'bit'
        @param cigar: a C{str} CIGAR string as present in the 6th field
            of a SAM file.
        @param MD: a C{str} MD string that may be present in the
            optional field of a SAM file.
        @raise ValueError: If category != 'seqID' or 'bit'
        @return: a C{float} representing sequence identity or a
            C{int} representing bitscore.
        """

        if category is not 'seqID' and category is not 'bit':
            raise ValueError('Category must be "seqID" or "bit"')

        if 'M' not in cigar:
            sepNumLet = [''.join(v) for k, v in groupby(cigar, str.isdigit)]
            els = {'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0,
                   'X': 0}
            index = 0
            for element in xrange(0, len(sepNumLet), 2):
                els[sepNumLet[index + 1]] += int(sepNumLet[index])
                index += 2
            if category is 'seqID':
                seqLen = els['I'] + els['S'] + els['='] + els['X']
                seqID = els['='] / float(seqLen)
                return seqID
            else:
                bit = (match * els['='] + mismatch * els['X']
                       + insert * els['I'] + delete * els['D'])
                return bit
        else:
            sepNumLet = [''.join(v) for k, v in groupby(cigar, str.isdigit)]
            els = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0,
                   '=': 0, 'X': 0}
            index = 0
            for element in xrange(0, len(sepNumLet), 2):
                els[sepNumLet[index + 1]] += int(sepNumLet[index])
                index += 2
            seqLen = els['M'] + els['I'] + els['S'] + els['='] + els['X']

            sepNumLetMD = [''.join(v) for k, v in groupby(MD, str.isdigit)]
            elsMD = {'matchMD': 0, 'misMD': 0, 'delMD': 0}
            for item in sepNumLetMD:
                if str.isdigit(item):
                    elsMD['matchMD'] += int(item)
                elif '^' in item:
                    elsMD['delMD'] += len(item) - 1
                else:
                    elsMD['misMD'] += len(item)
            seqLenMD = elsMD['matchMD'] + elsMD['misMD']

            assert (seqLen - seqLenMD == els['I']), "seq lengths from MD and cigar are not equal"
            assert (els['D'] == elsMD['delMD']), "no. of deletions from MD and cigar are not equal"

            if category is 'seqID':
                seqID = elsMD['matchMD'] / float(seqLen)
                return seqID
            else:
                bit = (match * elsMD['matchMD'] + mismatch * elsMD['misMD']
                       + insert * els['I'] + delete * elsMD['delMD'])
                return bit

    def _calcHsp(self, pos, sbjctemplateLen, query, flag, score, seqid):
        """
        Calculate the offsets of subject and query, so that they are 0-based,
        and usable in python.

        @param pos: a C{int}. The POS field as reported in a .SAM file.
            Corresponds to the '1-based leftmost mapping POSition'
        @param sbjctemplateLen: a C{int}. The length of the subject.
        @param queryLen: a C{int}. The length of the query.
        @param flag: C{int}, bitwise FLAG. For unpaired reads: 0 is single-end
            read mapped in the forward strand, 16 is single-end read mapped in
            the reverse strand, 4 is unmapped.
        """
        queryLen = len(query)
        subjectStart = 0
        subjectEnd = queryLen - 1  # not sure about this
        queryStart = pos - 1
        queryEnd = pos + queryLen - 1
        frame = [1, 1]
        if 'read reverse strand' in explain_sam_flags(flag):
            frame[0] = -1

        hsp = {"sbjct_end": subjectEnd,
               "expect": seqid,
               "sbjct": ' ',
               "sbjct_start": subjectStart,
               "query": query,
               "frame": frame,
               "query_end": queryEnd,
               "query_start": queryStart,
               "bits": score}

        return hsp

    def _convertSamAlignmentsToDict(self, line):
        """
        Pull (only) the fields we use out of the samFile and return them as a
        dict.

        @param line: an instance of a SAM file.

        Optional values after qual will vary between aligners and whether
        the read is aligned.
        """
        if not line.startswith('$') and not line.startswith('@'):
            line = line.strip().split()
            query = line[0]
            flag = int(line[1])
            refSeqName = line[2]
            pos = int(line[3])
            mapq = int(line[4])
            cigar = line[5]
            seq = line[9]
            qual = line[10]
            optional = line[11:]

            lengths = self.params['@SQ']  # This doesn't work

            if refSeqName != '*':
                if not cigar.isalnum() and '=' not in cigar:
                    raise ValueError("Cigar string is not alphanumeric - with "
                                     "the exception of =")

                subjTempLen = lengths[refSeqName]  # does this work? No

                record = {
                    "query": query,
                    "alignment": {
                        "length": subjTempLen,
                        "hsps": [],
                        "title": refSeqName,
                    }
                }

                optionalDict = {}
                for tag in optional:
                    key = tag.split(':')[0]
                    val = tag.split(':')[2]
                    optionalDict[key] = val

                if 'M' in cigar:
                    try:
                        MD = optionalDict['MD']
                        seqid = self._convertCigarMD('seqID', cigar, MD=MD)
                    except KeyError:
                        raise ValueError("No optional tag MD."
                                         "Run function findMD before"
                                         "proceeding")

                    try:
                        score = optionalDict['AS']
                    except KeyError:
                        score = self._convertCigarMD('bit', cigar, MD=MD)
                else:
                    seqid = self._convertCigarMD('seqID', cigar)
                    try:
                        score = optionalDict['AS']
                    except KeyError:
                        score = self._convertCigarMD('bit', cigar)

                record["alignment"]["hsps"].append(self._calcHsp(pos,
                                                                 subjTempLen,
                                                                 seq, flag,
                                                                 score,
                                                                 seqid))
            else:
                record = {
                    "query": line[0],
                }

            return record

    def _convertSamHeaderToDict(self, record):
        """
        Takes the lines of a SAM file beginning with '@' and returns a dict.
        For @SQ lines, the ref seq name is a key and the length the value.

        @param record: an instance of a SAM file which starts with '@'.
        @return: a C{dict}.
        """
        result = {'@HD': [], '@SQ': [], '@RG': [], '@PG': [], '@CO': []}
        for line in record:
            line = line.strip().split()
            if line[0] == '@SQ':
                tags = line[1:]
                tagDict = {}
                for tag in tags:
                    if 'SN' in tag:
                        key = tag.split(':')[1]
                    elif 'LN' in tag:
                        val = tag.split(':')[1]
                tagDict[key] = val
                result[line[0]].append(tagDict)
            elif line[0] in result:
                tags = line[1:]
                tagDict = {}
                for tag in tags:
                    key = tag.split(':')[0]
                    val = tag.split(':')[1]
                    tagDict[key] = val
                result[line[0]].append(tagDict)
        return result

    def records(self):
        """
        Yield SAM records. Set self.params from data in the SAM header.
        """
        if not self.filename.lower().endswith('.sam'):
            raise ValueError('A SAM file must be given')
        headerLines = []
        with open(self.filename) as fp:
            for record in fp:
                if record.startswith('@'):
                    headerLines.append(record)
                yield record

            self.params = self._convertSamHeaderToDict(headerLines)

    def saveAsJSON(self):
        """
        Write the records out as JSON. The first JSON object saved contains
        the SAM header tags.

        @param outFile: a C{str} name of the output file.
        """
        if not self.outFile.lower().endswith('.json'):
            raise ValueError('A JSON file must be given')
        with open(self.outFile, 'w') as ofp:
            first = True
            for rec in self.records():
                if rec.startswith('@') and first:
                    print >>ofp, dumps(self.params, separators=(',', ':'))
                    first = False
                elif not rec.startswith('@'):
                    print >>ofp, dumps(self._convertSamAlignmentsToDict(rec),
                                       separators=(',', ':'))


class SAMRecordsReader(object):
    """
    Provide a method that yields alignment records from a SAM file.
    Store, check, and make accessible the alignment parameters.

    @param filename: A C{str} filename containing alignment records.
    @param application: The C{str} name of the alignment program used.
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}.
    """

    def __init__(self, filename, applicationParams=None,
                 scoreClass='HigherIsBetterScore'):
        self._filename = filename
        self._scoreClass = scoreClass
        self._appParams = applicationParams
        if scoreClass is 'HigherIsBetterScore':
            self._hspClass = HSP
        # else:
            # self._hspClass = LSP

        self._open(filename)
        self.application = applicationParams['app'].lower()

    def _open(self, filename):
        """
        Open the input file. Set self._fp to point to it. Read the first
        line of parameters.

        @param filename: A C{str} filename containing alignment records.
        @raise ValueError: if the first line of the file isn't a valid SAM
            file.
        """
        if filename.lower().endswith('.sam'):
            self._fp = open(filename)
        else:
            raise ValueError('Invalid file type given: %s' % self._filename)

        line = self._fp.readline()
        if not line:
            raise ValueError('SAM file %r was empty.' % self._filename)

    def _convertCigarMD(self, cigar, MD=None, match=6,
                        mismatch=-3, insert=-1, delete=-1):
        """
        Calculates 'bitscore' from CIGAR and MD string.
        Bitscore scoring: =:6, X:-3, I:-1, D:-1, can be changed.

        @param cigar: a C{str} CIGAR string as present in the 6th field
            of a SAM file.
        @param MD: a C{str} MD string that may be present in the
            optional field of a SAM file.
        @raise ValueError: If category != 'seqID' or 'bit'
        @return: a C{int} representing bitscore.
        """
        if 'M' not in cigar:
            sepNumLet = [''.join(v) for k, v in groupby(cigar, str.isdigit)]
            els = {'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0,
                   'X': 0}
            index = 0
            for element in xrange(0, len(sepNumLet), 2):
                els[sepNumLet[index + 1]] += int(sepNumLet[index])
                index += 2
            bit = (match * els['='] + mismatch * els['X'] + insert * els['I']
                   + delete * els['D'])
            return bit

        else:
            sepNumLet = [''.join(v) for k, v in groupby(cigar, str.isdigit)]
            els = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0,
                   '=': 0, 'X': 0}
            index = 0
            for element in xrange(0, len(sepNumLet), 2):
                els[sepNumLet[index + 1]] += int(sepNumLet[index])
                index += 2
            seqLen = els['M'] + els['I'] + els['S'] + els['='] + els['X']

            sepNumLetMD = [''.join(v) for k, v in groupby(MD, str.isdigit)]
            elsMD = {'matchMD': 0, 'misMD': 0, 'delMD': 0}
            for item in sepNumLetMD:
                if str.isdigit(item):
                    elsMD['matchMD'] += int(item)
                elif '^' in item:
                    elsMD['delMD'] += len(item) - 1
                else:
                    elsMD['misMD'] += len(item)
            seqLenMD = elsMD['matchMD'] + elsMD['misMD']

            assert (seqLen - seqLenMD == els['I']), "seq lengths from MD and cigar are not equal"
            assert (els['D'] == elsMD['delMD']), "no. of deletions from MD and cigar are not equal"

            bit = (match * elsMD['matchMD'] + mismatch * elsMD['misMD']
                   + insert * els['I'] + delete * elsMD['delMD'])
            return bit

    def _lineToHSP(self, line):
        """
        Take a line of a SAM file and convert it to an HSP.

        @param line: Line of a SAM file.
        @raise ValueError: If the query id in the SAM entry does not
            match the id of the read.
        @return: A L{dark.alignment.Alignment.HSP} instance.
        """
        if not line.startswith('$') and not line.startswith('@'):
            line = line.strip().split()
            query = line[0]
            flag = int(line[1])
            refSeqName = line[2]
            pos = int(line[3])
            mapq = int(line[4])
            cigar = line[5]
            rnext = line[6]
            pnext = int(line[7])
            templateLen = abs(int(line[8]))
            seq = line[9]
            qual = line[10]
            optional = line[11:]

            if refSeqName != '*':
                if not cigar.isalnum() and '=' not in cigar:
                    raise ValueError("Cigar string is not alphanumeric - with"
                                     " the exception of =")

                optionalDict = {}
                for tag in optional:
                    key = tag.split(':')[0]
                    val = tag.split(':')[2]
                    optionalDict[key] = val
                try:
                    score = optionalDict['AS']
                # Need to make sure if 'AS' is used as a score for one read
                # then it's used for all of the others - can't have mix of
                # AS and score from _convertCigarMD
                except KeyError:
                    if 'M' in cigar:
                        try:
                            MD = optionalDict['MD']
                        except KeyError:
                            raise ValueError("No optional tag MD."
                                             "Run function findMD before "
                                             "proceeding")
                        else:
                            score = self._convertCigarMD(cigar, MD=MD)
                    else:
                        score = self._convertCigarMD(cigar)

                cigarSplit = [''.join(v) for k, v in groupby(cigar,
                                                             str.isdigit)]
                if cigarSplit[1] is 'S':
                    beginS = cigarSplit[0]
                elif cigarSplit[1] is 'H' and cigarSplit[3] is 'S':
                    beginS = cigarSplit[2]
                else:
                    beginS = 0

                if cigarSplit[-1] is 'S':
                    lastS = cigarSplit[-2]
                elif cigarSplit[-1] is 'H' and cigarSplit[-3] is 'S':
                    lastS = cigarSplit[-4]
                else:
                    lastS = 0

                cig = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0,
                       'P': 0, '=': 0, 'X': 0}
                index = 0
                for element in xrange(0, len(cigarSplit), 2):
                    cig[cigarSplit[index + 1]] += int(cigarSplit[index])
                    index += 2

                readStart = 0 + beginS
                readEnd = len(seq) - 1 - lastS
                subjStart = pos - 1
                subjEnd = (pos - 2 + len(seq) - cig['I'] - cig['S']
                           - cig['P'] + cig['N'] + cig['D'])

                hsp = HSP(score, readStart=readStart, readEnd=readEnd,
                          subjectStart=subjStart, subjectEnd=subjEnd)
                return hsp

    def readAlignments(self, reads):
        """
        Read lines of SAM file from self._filename, convert them to read
        alignments and yield them.

        @param reads: A generator yielding L{Read} instances, corresponding to
            the reads that were given to the aligner.
        @raise ValueError: If any of the lines in the file cannot be parsed.
        @return: A generator that yields C{dark.alignments.ReadAlignments}
            instances.
        """
        self._open(self._filename)  # Function to open the file.
        readIds = defaultdict(lambda : defaultdict(list))
        for line in self._fp:
            if not line.startswith('$') and not line.startswith('@'):
                details = line.strip().split()  # Gives list of strings
                query = details[0]
                refSeqName = details[2]
                if refSeqName != '*':
                    subjTempLen = self._appParams[refSeqName]  # Does this work?
                    alignment = Alignment(subjTempLen, refSeqName)
                    readIds[query][alignment].append(self._lineToHSP(line))

        for read in reads:
            try:
                alignments = readIds[read.id]
            except KeyError:
                raise ValueError(
                    'Read id %s found in passed reads but not in SAM file %r '
                    % (read.id, self._fp))
            else:
                yield ReadAlignments(read, alignments)

        self._fp.close()


def findMD(samFile, fastaFile):
    """
    Calculates the optional MD tag. If MD tag already present
    raises error if any are found to be different.
    Requires samtools.

    @param samFile: a C{str} of a SAM file.
    @param fastaFile: a C{str} of the FASTA file used for the
        SAM file.
    @return: samFile with the MD strings.
    """
    if not fastaFile.lower().endswith('.fasta'):
        raise ValueError('A FASTA file must be given')
    if not samFile.lower().endswith('.sam'):
        raise ValueError('A SAM file must be given')

    samFileNew = ''.join(samFile.split())[:-4] + '-fillmd.sam'
    sp.Popen(['samtools', 'fillmd', '-S', samFile, fastaFile, '>',
             samFileNew], stderr=sp.PIPE)
    return samFileNew

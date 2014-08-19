from collections import defaultdict

from dark.sam.conversion import SAMRecordsReader
from dark.score import HigherIsBetterScore
from dark.hsp import HSP
from dark.alignments import (
    Alignment, ReadsAlignments, ReadAlignments, ReadsAlignmentsParams)


class SamReadsAlignments(ReadsAlignments):
    """
    Hold information about a set of SAM records.

    @param reads: A L{dark.reads.Reads} instance providing the sequences that
        were given to the alignment application as queries. Note that the
        order of the reads *MUST* match the order of the records in the
        application output files.
    @param samFilename: Either a single C{str} filename or a C{list} of
        C{str} file names containing SAM output.
    @raises ValueError: if a file type is not recognized, if the number of
        reads does not match the number of records found in the BLAST result
        files, or if BLAST parameters in all files do not match.
    """

    def _convertSamHeaderToDict(self):
        """
        Takes the lines of a SAM file beginning with '@' and returns a dict.
        For @SQ lines, the ref seq name is a key and the length the value.

        @param record: a C{str} of samFilename.
        @return: a C{dict}.
        """
        headerLines = []
        with open(self.samFilename) as fp:
            for record in fp:
                if record.startswith('@'):
                    headerLines.append(record)
            result = {'@HD': [], '@SQ': [], '@RG': [], '@PG': [], '@CO': []}
            for line in headerLines:
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

    def __init__(self, reads, samFilename):
        self.samFilename = samFilename
        self.reads = reads
        # Prepare application parameters in order to initialize self.
        header = _convertSamHeaderToDict(self)
        app = header['@PG']['ID']

        applicationParams = ReadsAlignmentsParams('application'=app,
                                                  'samParams'=header)

        ReadsAlignments.__init__(self, reads, applicationParams,
                                 scoreClass=HigherIsBetterScore)

    def _getReader(self, filename, scoreClass):
        """
        Obtain a SAM record reader.

        @param filename: The C{str} file name holding the SAM records.
        @param scoreClass: A class to hold and compare scores (see scores.py).
        """
        if filename.endswith('.sam'):
            return SAMRecordsReader(samFilename, scoreClass)
        # Haven't written JSONRecordsReader for after SAM->JSON
        # elif filename.endswith('.json'):
            # return JSONRecordsReader(samFilename, scoreClass)
        else:
            raise ValueError(
                'Unknown SAM file suffix for file %r.' % filename)

    def iter(self):
        """
        Extract SAM records and yield C{ReadAlignments} instances.

        @return: A generator that yields C{ReadAlignments} instances.
        """

        reader = self._getReader(samFilename, self.scoreClass)

        for readAlignments in reader.readAlignments(reads):
            count += 1
            yield readAlignments

    # Is anything calling this function?
    def getSequence(self, title):
        """
        Obtain information about a sequence, given its title.

        @param title: A C{str} sequence title from a BLAST hit. Of the form
            'gi|63148399|gb|DQ011818.1| Description...'.
        @return: A C{SeqIO.read} instance.
        """
        # Look up the title in the database that was given to BLAST on the
        # command line.
        return ncbidb.getSequence(
            title, self.params.applicationParams['database'])

from dark.sam.conversion import SAMRecordsReader
from dark.score import HigherIsBetterScore
from dark.alignments import (ReadsAlignments, ReadsAlignmentsParams)


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
            if headerLines is None:
                raise ValueError('No header lines in %s' % self.samFilename)
            result = {}
            for line in headerLines:
                line = line.strip().split()
                tags = line[1:]
                if line[0] == '@SQ':
                    for tag in tags:
                        if 'SN' in tag:
                            key = tag.split(':')[1]
                        elif 'LN' in tag:
                            val = tag.split(':')[1]
                        result[key] = val
                elif line[0] == '@PG':
                    for tag in tags:
                        if 'ID' in tag:
                            result['app'] = tag.split(':')[1]
                        else:
                            key = tag.split(':')[0]
                            val = tag.split(':')[1]
                            result[key] = val
                else:
                    for tag in tags:
                        key = tag.split(':')[0]
                        val = tag.split(':')[1]
                        result[key] = val
            return result

    def __init__(self, reads, samFilename):
        self.samFilename = samFilename
        self.reads = reads
        # Prepare application parameters in order to initialize self.
        self.head = self._convertSamHeaderToDict()
        app = self.head['app']

        applicationParams = ReadsAlignmentsParams(app,
                                                  applicationParams=self.head)

        ReadsAlignments.__init__(self, reads, applicationParams,
                                 scoreClass=HigherIsBetterScore)

    def _getReader(self, filename, applicationParams, scoreClass):
        """
        Obtain a SAM record reader.

        @param filename: The C{str} file name holding the SAM records.
        @param scoreClass: A class to hold and compare scores (see scores.py).
        """
        if filename.endswith('.sam'):
            return SAMRecordsReader(self.samFilename, self.head,
                                    self.scoreClass)
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
        count = 0
        reads = iter(self.reads)
        reader = self._getReader(self.samFilename, self.head,
                                 self.scoreClass)

        for readAlignments in reader.readAlignments(reads):
            count += 1
            yield readAlignments

        # Make sure all reads were used.
        try:
            read = reads.next()
        except StopIteration:
            pass
        else:
            raise ValueError(
                'Reads iterator contained more reads than the number of '
                'BLAST records found (%d). First unknown read id is %r.' %
                (count, read.id))

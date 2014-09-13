from dark.sam.conversion import SAMRecordsReader
from dark.score import HigherIsBetterScore
from dark.alignments import ReadsAlignments, ReadsAlignmentsParams
from dark.sam.hacks import checkSAMfile


class SamReadsAlignments(ReadsAlignments):
    """
    Hold information about a set of SAM records.
    TODO If there are query names that come in pairs then they're paired-end
        reads. Messes up the principle of the new code - the paired
        end reads have the same identifier but are different sequences,
        therefore counting them as two HSPs of one read alignment is not
        accurate.

    @param reads: A L{dark.reads.Reads} instance providing the sequences that
        were given to the alignment application as queries. Note that the
        order of the reads *MUST* match the order of the records in the
        application output files.
    @param samFilename: Either a single C{str} filename or a C{list} of
        C{str} file names containing SAM output.
    @raises ValueError: if a file type is not recognized, if the number of
        reads does not match the number of records found in the SAM
        files.
    """
    def __init__(self, reads, samFilename):
        if checkSAMfile(samFilename):
            self.samFilename = samFilename
        self.reads = reads
        # Prepare application parameters in order to initialize self.
        self.head = self._convertSamHeaderToDict()
        app = self.head['application']

        applicationParams = ReadsAlignmentsParams(app,
                                                  applicationParams=self.head)

        ReadsAlignments.__init__(self, reads, applicationParams,
                                 scoreClass=HigherIsBetterScore)

    def _convertSamHeaderToDict(self):
        """
        Takes the lines of a SAM file beginning with '@' and returns a dict.
        For @SQ lines, the ref seq name is a key and the length the value.

        @param record: a C{str} of samFilename.
        @return: a C{dict} of the parameters found in the header of a
            SAMRecordsReader file.
        """
        headerLines = []
        with open(self.samFilename) as fp:
            for record in fp:
                if record[0] == '@':
                    headerLines.append(record)
                else:
                    break
        if not headerLines:
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
                        val = int(tag.split(':')[1])
                    result[key] = val
            elif line[0] == '@PG':
                for tag in tags:
                    if 'ID' in tag:
                        result['application'] = tag.split(':')[1]
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

    def _getReader(self, filename, applicationParams, scoreClass):
        """
        Obtain a SAM record reader.

        @param filename: The C{str} file name holding the SAM records.
        @param scoreClass: A class to hold and compare scores (see scores.py).
        """
        return SAMRecordsReader(self.samFilename, self.head,
                                self.scoreClass)

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
                'SAM records found (%d). First unknown read id is %r.' %
                (count, read.id))

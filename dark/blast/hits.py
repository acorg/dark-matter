import copy
import string

from dark.hits import Hit, Hits
from dark.blast.conversion import JSONRecordsReader
from dark.blast.params import checkCompatibleParams


def numericallySortFilenames(names):
    """
    Sort (ascending) a list of file names by their numerical prefixes.

    @param: A C{list} of file names, each of which starts with a string of
        digits.

    @return: The sorted C{list}.
    """

    def extractNumericPrefix(name):
        """
        Find any numeric prefix at the start of C{name} and return it as an
        C{int}.

        @param: A C{str} file name, possibly starting with some digits.
        @return: The C{int} number at the start of the name, else 0 if there
            are no leading digits.
        """
        count = 0
        for ch in name:
            if ch in string.digits:
                count += 1
            else:
                break
        return 0 if count == 0 else int(name[0:count])

    return sorted(names, key=lambda name: extractNumericPrefix(name))


class BlastHits(Hits):
    """
    Hold information about a set of BLAST records.

    @param reads: A L{dark.reads.Reads} instance providing the sequences that
        were given to BLAST as queries. Note that the order of the reads
        *MUST* match the order of the records in the BLAST output files.
    @param blastFilenames: Either a single C{str} filename or a C{list} of
        C{str} file names containing BLAST output. Files can either be XML
        (-outfmt 5) BLAST output file or our smaller (possibly bzip2
        compressed) converted JSON equivalent produced by
        C{bin/convert-blast-xml-to-json.py} from a BLAST XML file.
    @param scoreType: A C{str}, either 'bits' or 'e values'.
    @param sortBlastFilenames: A C{bool}. If C{True}, C{blastFilenames} will be
        sorted by numeric prefix (using L{numericallySortFilenames}) before
        being read. This can be used to conveniently sort the files produced
        by our HTCondor jobs.
    @raises ValueError: If an invalid C{scoreType} is passed, if a file type
        is not recognized, if the number of reads does not match the number
        of records found in the BLAST result files, or if BLAST parameters
        in all files do not match.
    """

    def __init__(self, reads, blastFilenames, scoreType='bits',
                 sortBlastFilenames=True):
        if type(blastFilenames) == str:
            blastFilenames = [blastFilenames]
        if sortBlastFilenames:
            self.blastFilenames = numericallySortFilenames(blastFilenames)
        else:
            self.blastFilenames = blastFilenames

        if scoreType not in ('bits', 'e values'):
            raise ValueError("ScoreType parameter %r is invalid. Valid "
                             "options are 'bits' and 'e values'." % scoreType)
        self.scoreType = scoreType

        # Read and copy the BLAST parameters and initialize self.
        self._reader = self._getReader(self.blastFilenames[0])
        application = self._reader.params['application'].lower()
        Hits.__init__(self, reads, application,
                      copy.deepcopy(self._reader.params))

    def _getReader(self, filename):
        if filename.endswith('.json') or filename.endswith('.json.bz2'):
            return JSONRecordsReader(filename)
        else:
            raise ValueError(
                'Unknown BLAST record file suffix for file %r.' % filename)

    def __iter__(self):
        """
        Extract BLAST records and yield hits.

        For each file except the first, check that the BLAST parameters are
        compatible with those found (above, in __init__) in the first file.
        """
        # Note that self._reader is already initialized (in __init__) for
        # the first input file. This is less clean than it could be, but it
        # makes testing easier, since open() is then only called once for
        # each input file.

        count = 0
        reader = self._reader
        reads = iter(self.reads)
        first = True

        for blastFilename in self.blastFilenames:
            if first:
                first = False
            else:
                reader = self._getReader(blastFilename)
                differences = checkCompatibleParams(self.params, reader.params)
                if differences:
                    raise ValueError(
                        'Incompatible BLAST parameters found. The parameters '
                        'in %s differ from those originally found in %s. %s' %
                        (blastFilename, self.blastFilenames[0], differences))

            for read, alignments in reader.records(
                    reads, self.scoreType, self.application):
                count += 1
                yield Hit(read, alignments)

            reader.close()

        # Make sure all reads were used.
        try:
            read = reads.next()
        except StopIteration:
            pass
        else:
            raise ValueError(
                'Reads iterator contained more reads than the number of BLAST '
                'records found (%d). First unknown read id is %r.' %
                (count, read.id))

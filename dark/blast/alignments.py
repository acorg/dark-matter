import copy
import string

from dark.score import HigherIsBetterScore
from dark.alignments import ReadsAlignments, ReadsAlignmentsParams
from dark.blast.conversion import JSONRecordsReader
from dark.blast.params import checkCompatibleParams
from dark import ncbidb


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


class BlastReadsAlignments(ReadsAlignments):
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
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}, for comparing bit scores. If you
        are using e-values, pass LowerIsBetterScore instead.
    @param sortBlastFilenames: A C{bool}. If C{True}, C{blastFilenames} will be
        sorted by numeric prefix (using L{numericallySortFilenames}) before
        being read. This can be used to conveniently sort the files produced
        by our HTCondor jobs.
    @param randomizeZeroEValues: If C{True}, e-values that are zero will be set
        to a random (very good) value.
    @raises ValueError: if a file type is not recognized, if the number of
        reads does not match the number of records found in the BLAST result
        files, or if BLAST parameters in all files do not match.
    """

    def __init__(self, reads, blastFilenames, scoreClass=HigherIsBetterScore,
                 sortBlastFilenames=True, randomizeZeroEValues=True):
        if type(blastFilenames) == str:
            blastFilenames = [blastFilenames]
        if sortBlastFilenames:
            self.blastFilenames = numericallySortFilenames(blastFilenames)
        else:
            self.blastFilenames = blastFilenames
        self.randomizeZeroEValues = randomizeZeroEValues

        # Prepare application parameters in order to initialize self.
        self._reader = self._getReader(self.blastFilenames[0], scoreClass)
        application = self._reader.application
        blastParams = copy.deepcopy(self._reader.params)
        subjectIsNucleotides = application != 'blastx'
        scoreTitle = ('Bit score' if scoreClass is HigherIsBetterScore
                      else '$- log_{10}(e)')

        applicationParams = ReadsAlignmentsParams(
            application, blastParams,
            subjectIsNucleotides=subjectIsNucleotides, scoreTitle=scoreTitle)

        ReadsAlignments.__init__(self, reads, applicationParams,
                                 scoreClass=scoreClass)

    def _getReader(self, filename, scoreClass):
        """
        Obtain a JSON record reader for BLAST records.

        @param filename: The C{str} file name holding the JSON.
        @param scoreClass: A class to hold and compare scores (see scores.py).
        """
        if filename.endswith('.json') or filename.endswith('.json.bz2'):
            return JSONRecordsReader(filename, scoreClass)
        else:
            raise ValueError(
                'Unknown BLAST record file suffix for file %r.' % filename)

    def __iter__(self):
        """
        Extract BLAST records and yield C{ReadAlignments} instances.

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
                # No need to check params in the first file. We already read
                # them in and stored them in __init__.
                first = False
            else:
                reader = self._getReader(blastFilename, self.scoreClass)
                differences = checkCompatibleParams(
                    self.params.applicationParams, reader.params)
                if differences:
                    raise ValueError(
                        'Incompatible BLAST parameters found. The parameters '
                        'in %s differ from those originally found in %s. %s' %
                        (blastFilename, self.blastFilenames[0], differences))

            for readAlignments in reader.readAlignments(reads):
                count += 1
                yield readAlignments

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

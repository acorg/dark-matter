from random import uniform
from math import log10
import copy

from dark.alignments import (
    ReadsAlignments, ReadAlignments, ReadsAlignmentsParams)
from dark.diamond.conversion import JSONRecordsReader
from dark.fasta import FastaReads, SqliteIndex
from dark.reads import AAReadWithX
from dark.score import HigherIsBetterScore
from dark.utils import numericallySortFilenames

ZERO_EVALUE_UPPER_RANDOM_INCREMENT = 150


class DiamondReadsAlignments(ReadsAlignments):
    """
    Hold information about a set of DIAMOND records.

    @param reads: A L{dark.reads.Reads} instance providing the sequences that
        were given to DIAMOND as queries. Note that the order of the reads
        *MUST* match the order of the records in the DIAMOND output files.
    @param filenames: Either a single C{str} filename or a C{list} of C{str}
        file names containing our (possibly bzip2 compressed) per-line JSON
        produced by C{bin/convert-diamond-to-json.py} from DIAMOND XML output.
    @param filenames: Either a single C{str} filename or a C{list} of C{str}
        file names containing our (possibly bzip2 compressed) per-line JSON
        produced by C{bin/convert-diamond-to-json.py} from DIAMOND XML output.
    @param databaseFilename: A C{str} holding the name of the FASTA file used
        to make the DIAMOND database. Cannot be used with
        C{sqliteDatabaseFilename}.
    @param databaseDirectory: The directory where the FASTA file
        used to make the DIAMOND database can be found. This argument is only
        useful when sqliteDatabaseFilename is specified.
    @param sqliteDatabaseFilename: A C{str} holding the name of the sqlite3
        database file made from the FASTA used to make the DIAMOND database
        (Use ../../bin/make-fasta-database.py to construct such a database).
        Cannot be used with C{databaseFilename}.
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}, for comparing bit scores. If you
        are using e-values, pass LowerIsBetterScore instead.
    @param sortFilenames: A C{bool}. If C{True}, C{filenames} will be
        sorted by numeric prefix (using L{numericallySortFilenames}) before
        being read. This can be used to conveniently sort the files produced
        by our HTCondor jobs.
    @param randomizeZeroEValues: If C{True}, e-values that are zero will be set
        to a random (very good) value.
    @raises ValueError: if a file type is not recognized, or if the number of
        reads does not match the number of records found in the DIAMOND result
        files, or if neither (or both) of databaseFilename and
        sqliteDatabaseFilename are given.
    """

    def __init__(self, reads, filenames, databaseFilename=None,
                 databaseDirectory=None, sqliteDatabaseFilename=None,
                 scoreClass=HigherIsBetterScore, sortFilenames=True,
                 randomizeZeroEValues=True):
        if databaseFilename is None and sqliteDatabaseFilename is None:
            raise ValueError(
                'Either databaseFilename or sqliteDatabaseFilename must be '
                'provided to %s' % self.__class__.__name__)
        elif not (databaseFilename is None or sqliteDatabaseFilename is None):
            raise ValueError(
                'databaseFilename and sqliteDatabaseFilename cannot both be '
                'provided to %s' % self.__class__.__name__)
        if type(filenames) == str:
            filenames = [filenames]
        if sortFilenames:
            self.filenames = numericallySortFilenames(filenames)
        else:
            self.filenames = filenames

        self._databaseFilename = databaseFilename
        self._sqliteDatabaseFilename = sqliteDatabaseFilename
        self._databaseDirectory = databaseDirectory
        self._subjectTitleToSubject = None
        self.randomizeZeroEValues = randomizeZeroEValues

        # Prepare diamondTask parameters in order to initialize self.
        self._reader = self._getReader(self.filenames[0], scoreClass)
        diamondTask = self._reader.diamondTask
        diamondParams = copy.deepcopy(self._reader.params)
        scoreTitle = ('Bit score' if scoreClass is HigherIsBetterScore
                      else '$- log_{10}(e)$')

        diamondTaskParams = ReadsAlignmentsParams(
            diamondTask, diamondParams,
            subjectIsNucleotides=False,  # DIAMOND dbs are always protein.
            scoreTitle=scoreTitle)

        ReadsAlignments.__init__(self, reads, diamondTaskParams,
                                 scoreClass=scoreClass)

    def _getReader(self, filename, scoreClass):
        """
        Obtain a JSON record reader for DIAMOND records.

        @param filename: The C{str} file name holding the JSON.
        @param scoreClass: A class to hold and compare scores (see scores.py).
        """
        if filename.endswith('.json') or filename.endswith('.json.bz2'):
            return JSONRecordsReader(filename, scoreClass)
        else:
            raise ValueError(
                'Unknown DIAMOND record file suffix for file %r.' % filename)

    def iter(self):
        """
        Extract DIAMOND records and yield C{ReadAlignments} instances.

        @return: A generator that yields C{ReadAlignments} instances.
        """
        # Note that self._reader is already initialized (in __init__) for
        # the first input file. This is less clean than it could be, but it
        # makes testing easier, since open() is then only called once for
        # each input file.

        reads = iter(self.reads)
        first = True

        for filename in self.filenames:
            if first:
                # The first file has already been opened, in __init__.
                first = False
                reader = self._reader
            else:
                reader = self._getReader(filename, self.scoreClass)

            for readAlignments in reader.readAlignments(reads):
                yield readAlignments

        # Any remaining query reads must have had no subject matches.
        for read in reads:
            yield ReadAlignments(read, [])

    def getSubjectSequence(self, title):
        """
        Obtain information about a subject sequence given its title.

        @param title: A C{str} sequence title from a DIAMOND hit.
        @raise KeyError: If the C{title} is not present in the DIAMOND
            database.
        @return: An C{AAReadWithX} instance.
        """
        if self._subjectTitleToSubject is None:
            if self._databaseFilename is None:
                # An Sqlite3 database is used to look up subjects.
                self._subjectTitleToSubject = SqliteIndex(
                    self._sqliteDatabaseFilename,
                    fastaDirectory=self._databaseDirectory,
                    readClass=AAReadWithX)
            else:
                # Build a dict to look up subjects.
                titles = {}
                for read in FastaReads(self._databaseFilename,
                                       readClass=AAReadWithX):
                    titles[read.id] = read
                self._subjectTitleToSubject = titles

        return self._subjectTitleToSubject[title]

    def adjustHspsForPlotting(self, titleAlignments):
        """
        Our HSPs are about to be plotted. If we are using e-values, these need
        to be adjusted.

        @param titleAlignments: An instance of L{TitleAlignment}.
        """
        # If we're using bit scores, there's nothing to do.
        if self.scoreClass is HigherIsBetterScore:
            return

        # Convert all e-values to high positive values, and keep track of the
        # maximum converted value.
        maxConvertedEValue = None
        zeroHsps = []

        # Note: don't call self.hsps() here because that will read them
        # from disk again, which is not what's wanted.
        for hsp in titleAlignments.hsps():
            if hsp.score.score == 0.0:
                zeroHsps.append(hsp)
            else:
                convertedEValue = -1.0 * log10(hsp.score.score)
                hsp.score.score = convertedEValue
                if (maxConvertedEValue is None or
                        convertedEValue > maxConvertedEValue):
                    maxConvertedEValue = convertedEValue

        if zeroHsps:
            # Save values so that we can use them in self.adjustPlot
            self._maxConvertedEValue = maxConvertedEValue
            self._zeroEValueFound = True

            # Adjust all zero e-value HSPs to have numerically high values.
            if self.randomizeZeroEValues:
                for hsp in zeroHsps:
                    hsp.score.score = (maxConvertedEValue + 2 + uniform(
                        0, ZERO_EVALUE_UPPER_RANDOM_INCREMENT))
            else:
                for count, hsp in enumerate(zeroHsps, start=1):
                    hsp.score.score = maxConvertedEValue + count
        else:
            self._zeroEValueFound = False

    def adjustPlot(self, readsAx):
        """
        Add a horizontal line to the plotted reads if we're plotting e-values
        and a zero e-value was found.

        @param readsAx: A Matplotlib sub-plot instance, as returned by
            matplotlib.pyplot.subplot.
        """
        # If we're using bit scores, there's nothing to do.
        if self.scoreClass is HigherIsBetterScore:
            return

        if self._zeroEValueFound:
            readsAx.axhline(y=self._maxConvertedEValue + 0.5, color='#cccccc',
                            linewidth=0.5)

from random import uniform
from math import log10
import copy
from typing import Generator

from dark.score import HigherIsBetterScore
from dark.alignments import ReadsAlignments, ReadsAlignmentsParams, ReadAlignments
from dark.blast.conversion import JSONRecordsReader
from dark.blast.params import checkCompatibleParams
from dark.fasta import FastaReads, SqliteIndex
from dark.reads import AARead, DNARead
from dark.utils import numericallySortFilenames

ZERO_EVALUE_UPPER_RANDOM_INCREMENT = 150


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
    @param databaseFilename: A C{str} holding the name of the FASTA file used
        to make the BLAST database. Cannot be used with
        C{sqliteDatabaseFilename}.
    @param databaseDirectory: The directory where the FASTA file
        used to make the BLAST database can be found. This argument is only
        useful when sqliteDatabaseFilename is specified.
    @param sqliteDatabaseFilename: A C{str} holding the name of the sqlite3
        database file made from the FASTA used to make the BLAST database
        (Use ../../bin/make-fasta-database.py to construct such a database).
        Cannot be used with C{databaseFilename}.
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

    def __init__(
        self,
        reads,
        blastFilenames,
        databaseFilename=None,
        databaseDirectory=None,
        sqliteDatabaseFilename=None,
        scoreClass=HigherIsBetterScore,
        sortBlastFilenames=True,
        randomizeZeroEValues=True,
    ):
        if type(blastFilenames) == str:
            blastFilenames = [blastFilenames]
        if sortBlastFilenames:
            self.blastFilenames = numericallySortFilenames(blastFilenames)
        else:
            self.blastFilenames = blastFilenames
        self._databaseFilename = databaseFilename
        self._sqliteDatabaseFilename = sqliteDatabaseFilename
        self._databaseDirectory = databaseDirectory
        self._subjectTitleToSubject = None
        self.randomizeZeroEValues = randomizeZeroEValues

        # Prepare application parameters in order to initialize self.
        self._reader = self._getReader(self.blastFilenames[0], scoreClass)
        application = self._reader.application
        blastParams = copy.deepcopy(self._reader.params)
        subjectIsNucleotides = application != "blastx"
        scoreTitle = (
            "Bit score" if scoreClass is HigherIsBetterScore else "$- log_{10}(e)$"
        )

        applicationParams = ReadsAlignmentsParams(
            application,
            blastParams,
            subjectIsNucleotides=subjectIsNucleotides,
            scoreTitle=scoreTitle,
        )

        ReadsAlignments.__init__(self, reads, applicationParams, scoreClass=scoreClass)

    def _getReader(self, filename, scoreClass):
        """
        Obtain a JSON record reader for BLAST records.

        @param filename: The C{str} file name holding the JSON.
        @param scoreClass: A class to hold and compare scores (see scores.py).
        """
        if filename.endswith(".json") or filename.endswith(".json.bz2"):
            return JSONRecordsReader(filename, scoreClass)
        else:
            raise ValueError("Unknown BLAST record file suffix for file %r." % filename)

    def iter(self) -> Generator[ReadAlignments, None, None]:
        """
        Extract BLAST records and yield C{ReadAlignments} instances.

        For each file except the first, check that the BLAST parameters are
        compatible with those found (above, in __init__) in the first file.

        @return: A generator that yields C{ReadAlignments} instances.
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
                    self.params.applicationParams, reader.params
                )
                if differences:
                    raise ValueError(
                        "Incompatible BLAST parameters found. The parameters "
                        "in %s differ from those originally found in %s. %s"
                        % (blastFilename, self.blastFilenames[0], differences)
                    )

            for readAlignments in reader.readAlignments(reads):
                count += 1
                yield readAlignments

        # Make sure all reads were used.
        try:
            read = next(reads)
        except StopIteration:
            pass
        else:
            raise ValueError(
                "Reads iterator contained more reads than the number of BLAST "
                "records found (%d). First unknown read id is %r." % (count, read.id)
            )

    def getSubjectSequence(self, title):
        """
        Obtain information about a subject sequence given its title.

        This information is cached in self._subjectTitleToSubject. It can
        be obtained from either a) an sqlite database (given via the
        sqliteDatabaseFilename argument to __init__), b) the FASTA that was
        originally given to BLAST (via the databaseFilename argument), or
        c) from the BLAST database using blastdbcmd (which can be
        unreliable - occasionally failing to find subjects that are in its
        database).

        @param title: A C{str} sequence title from a BLAST hit. Of the form
            'gi|63148399|gb|DQ011818.1| Description...'.
        @return: An C{AARead} or C{DNARead} instance, depending on the type of
            BLAST database in use.

        """
        if self.params.application in {"blastp", "blastx"}:
            readClass = AARead
        else:
            readClass = DNARead

        if self._subjectTitleToSubject is None:
            if self._databaseFilename is None:
                if self._sqliteDatabaseFilename is None:
                    # Fall back to blastdbcmd.  ncbidb has to be imported
                    # as below so ncbidb.getSequence can be patched by our
                    # test suite.
                    from dark import ncbidb

                    seq = ncbidb.getSequence(
                        title, self.params.applicationParams["database"]
                    )
                    return readClass(seq.description, str(seq.seq))
                else:
                    # An Sqlite3 database is used to look up subjects.
                    self._subjectTitleToSubject = SqliteIndex(
                        self._sqliteDatabaseFilename,
                        fastaDirectory=self._databaseDirectory,
                        readClass=readClass,
                    )
            else:
                # Build an in-memory dict to look up subjects. This only
                # works for small databases, obviously.
                titles = {}
                for read in FastaReads(self._databaseFilename, readClass=readClass):
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
                if maxConvertedEValue is None or convertedEValue > maxConvertedEValue:
                    maxConvertedEValue = convertedEValue

        if zeroHsps:
            # Save values so that we can use them in self.adjustPlot
            self._maxConvertedEValue = maxConvertedEValue
            self._zeroEValueFound = True

            # Adjust all zero e-value HSPs to have numerically high values.
            if self.randomizeZeroEValues:
                for hsp in zeroHsps:
                    hsp.score.score = (
                        maxConvertedEValue
                        + 2
                        + uniform(0, ZERO_EVALUE_UPPER_RANDOM_INCREMENT)
                    )
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
            readsAx.axhline(
                y=self._maxConvertedEValue + 0.5, color="#cccccc", linewidth=0.5
            )

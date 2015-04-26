from random import uniform
from math import log10
import copy

from dark.score import HigherIsBetterScore
from dark.alignments import ReadsAlignments, ReadsAlignmentsParams
from dark.blast.conversion import JSONRecordsReader
from dark.blast.params import checkCompatibleParams
from dark import ncbidb
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
                      else '$- log_{10}(e)$')

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

    def iter(self):
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
                    self.params.applicationParams, reader.params)
                if differences:
                    raise ValueError(
                        'Incompatible BLAST parameters found. The parameters '
                        'in %s differ from those originally found in %s. %s' %
                        (blastFilename, self.blastFilenames[0], differences))

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
                'Reads iterator contained more reads than the number of BLAST '
                'records found (%d). First unknown read id is %r.' %
                (count, read.id))

    def getSubjectSequence(self, title):
        """
        Obtain information about a subject sequence, given its title.

        @param title: A C{str} sequence title from a BLAST hit. Of the form
            'gi|63148399|gb|DQ011818.1| Description...'.
        @return: A C{SeqIO.read} instance.
        """
        # Look up the title in the database that was given to BLAST on the
        # command line.
        return ncbidb.getSequence(
            title, self.params.applicationParams['database'])

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

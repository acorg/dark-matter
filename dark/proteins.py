from __future__ import division, print_function

import re
from os.path import dirname, join
from operator import itemgetter
from six.moves.urllib.parse import quote
import numpy as np
from textwrap import fill
from collections import Counter

try:
    import matplotlib.pyplot as plt
except ImportError:
    import platform
    if platform.python_implementation() == 'PyPy':
        raise NotImplementedError(
            'matplotlib is not supported under pypy')
    else:
        raise

from dark.dimension import dimensionalIterator
from dark.fasta import FastaReads
from dark.fastq import FastqReads
from dark.html import NCBISequenceLinkURL
from dark.reads import Reads

# The following regex is deliberately greedy (using .*) to consume the
# whole protein title before backtracking to find the last [virus title]
# section. That way, it will match just the last [virus title] in a
# protein. This avoids situations in which two [...] delimited substrings
# are present in a protein name (in which case we just want the last).
# E.g., the following is a complete protein name:
#
#   gi|19919894|ref|NP_612577.1| Enzymatic polyprotein [Contains: Aspartic
#   protease; Endonuclease; Reverse transcriptase] [Carnation etched ring
#   virus]
#
# Unfortunately the regex doesn't find the virus title when the protein
# title has nested [...] sections, as in this example:
#
#   gi|224808893|ref|YP_002643049.1| replication-associated protein [Tomato
#   leaf curl Nigeria virus-[Nigeria:2006]]
#
# I decided not to worry about nested [...] sections (there are only 2
# instances that I know of).
_VIRUS_RE = re.compile('^(.*)\[([^\]]+)\]$')

# The virus title assigned to proteins whose title strings cannot be parsed
# for a virus title (see previous comment).  Do not use '<', '>' or any
# other HTML special chars in the following.
_NO_VIRUS_TITLE = '[no virus name found in protein title]'


def splitTitles(titles):
    """
    Split a title like "Protein name [virus name]" into two pieces using
    the final square brackets to delimit the virus name.

    @param titles: A C{str} protein and virus title.
    @return: A 2-C{tuple} giving the C{str} protein title and C{str} virus
        title. If C{titles} cannot be split on square brackets, it is
        returned as the first tuple element, followed by _NO_VIRUS_TITLE.
    """
    match = _VIRUS_RE.match(titles)
    if match:
        proteinTitle = match.group(1).strip()
        virusTitle = match.group(2).strip()
    else:
        proteinTitle = titles
        virusTitle = _NO_VIRUS_TITLE

    return proteinTitle, virusTitle


def getVirusProteinCounts(filename):
    """
    Get the number of proteins for each virus in C{filename}.

    @param filename: Either C{None} or a C{str} FASTA file name. If C{None}
        an empty C{Counter} is returned. If a FASTA file name is given, its
        sequence ids should have the format used in the NCBI viral protein
        file, in which the protein name is followed by the virus name in
        square brackets.
    @return: A C{Counter} keyed by C{str} virus name, whose values are C{int}s
        with the count of the number of proteins for the virus.
    """
    result = Counter()
    if filename:
        for protein in FastaReads(filename):
            _, virusTitle = splitTitles(protein.id)
            if virusTitle != _NO_VIRUS_TITLE:
                result[virusTitle] += 1

    return result


class VirusSampleFiles(object):
    """
    Maintain a cache of virus/sample FASTA/FASTQ file names, creating
    de-duplicated FASTA/FASTQ files (from reads for all proteins of a virus
    that a sample has), on demand.

    @param proteinGrouper: An instance of C{ProteinGrouper}.
    @param format_: A C{str}, either 'fasta' or 'fastq' indicating the format
        of the files containing the reads matching proteins.
    @raise ValueError: If C{format_} is unknown.
    """
    def __init__(self, proteinGrouper, format_='fasta'):
        self._proteinGrouper = proteinGrouper
        if format_ in ('fasta', 'fastq'):
            self._format = format_
            self._readsClass = FastaReads if format_ == 'fasta' else FastqReads
        else:
            raise ValueError("format_ must be either 'fasta' or 'fastq'.")
        self._viruses = {}
        self._samples = {}
        self._readsFilenames = {}

    def add(self, virusTitle, sampleName):
        """
        Add a virus title, sample name combination and get its FASTA/FASTQ file
        name and unique read count. Write the FASTA/FASTQ file if it does not
        already exist. Save the unique read count into
        C{self._proteinGrouper}.

        @param virusTitle: A C{str} virus title.
        @param sampleName: A C{str} sample name.
        @return: A C{str} giving the FASTA/FASTQ file name holding all the
            reads (without duplicates) from the sample that matched the
            proteins in the given virus.
        """
        virusIndex = self._viruses.setdefault(virusTitle, len(self._viruses))
        sampleIndex = self._samples.setdefault(sampleName, len(self._samples))

        try:
            return self._readsFilenames[(virusIndex, sampleIndex)]
        except KeyError:
            reads = Reads()
            for proteinMatch in self._proteinGrouper.virusTitles[
                    virusTitle][sampleName]['proteins'].values():
                for read in self._readsClass(proteinMatch['readsFilename']):
                    reads.add(read)
            saveFilename = join(
                proteinMatch['outDir'],
                'virus-%d-sample-%d.%s' % (virusIndex, sampleIndex,
                                           self._format))
            reads.filter(removeDuplicates=True)
            nReads = reads.save(saveFilename, format_=self._format)
            # Save the unique read count into self._proteinGrouper
            self._proteinGrouper.virusTitles[
                virusTitle][sampleName]['uniqueReadCount'] = nReads
            self._readsFilenames[(virusIndex, sampleIndex)] = saveFilename
            return saveFilename

    def lookup(self, virusTitle, sampleName):
        """
        Look up a virus title, sample name combination and get its FASTA/FASTQ
        file name and unique read count.

        This method should be used instead of C{add} in situations where
        you want an exception to be raised if a virus/sample combination has
        not already been passed to C{add}.

        @param virusTitle: A C{str} virus title.
        @param sampleName: A C{str} sample name.
        @raise KeyError: If the virus title or sample name have not been seen,
            either individually or in combination.
        @return: A (C{str}, C{int}) tuple retrieved from self._readsFilenames
        """
        virusIndex = self._viruses[virusTitle]
        sampleIndex = self._samples[sampleName]
        return self._readsFilenames[(virusIndex, sampleIndex)]


class ProteinGrouper(object):
    """
    Group matched proteins by the virus they come from.

    @param assetDir: The C{str} directory name where
        C{noninteractive-alignment-panel.py} put its HTML, blue plot and
        alignment panel images, and FASTA or FASTQ files. This must be relative
        to the filenames that will later be passed to C{addFile}.
    @param sampleNameRegex: A C{str} regular expression that can be used to
        extract a short sample name from full file names subsequently passed
        to C{self.addFile}. The regular expression must have a matching group
        (delimited by parentheses) to capture the part of the file name that
        should be used as the sample name.
    @param format_: A C{str}, either 'fasta' or 'fastq' indicating the format
        of the files containing the reads matching proteins.
    @param proteinFastaFilename: If not C{None}, a C{str} filename giving the
        name of the FASTA file with the protein AA sequences with their
        associated viruses in square brackets. This is the format used by NCBI
        for the viral protein files. If given, the contents of this file will
        be used to determine how many proteins each matched virus has.
    @raise ValueError: If C{format_} is unknown.
    """

    VIRALZONE = 'http://viralzone.expasy.org/cgi-bin/viralzone/search?query='

    def __init__(self, assetDir='out', sampleNameRegex=None, format_='fasta',
                 proteinFastaFilename=None):
        self._assetDir = assetDir
        self._sampleNameRegex = (re.compile(sampleNameRegex) if sampleNameRegex
                                 else None)
        if format_ in ('fasta', 'fastq'):
            self._format = format_
        else:
            raise ValueError("format_ must be either 'fasta' or 'fastq'.")

        self._virusProteinCount = getVirusProteinCounts(proteinFastaFilename)

        # virusTitles will be a dict of dicts of dicts. The first two keys
        # will be a virus title and a sample name. The final dict will
        # contain 'proteins' (a list of dicts) and 'uniqueReadCount' (an int).
        self.virusTitles = {}
        # sampleNames is keyed by sample name and will have values that hold
        # the sample's alignment panel index.html file.
        self.sampleNames = {}
        self.virusSampleFiles = VirusSampleFiles(self, format_=format_)

    def _title(self):
        """
        Create a title summarizing the viruses and samples.

        @return: A C{str} title.
        """
        return (
            'Overall, proteins from %d virus%s were found in %d sample%s.' %
            (len(self.virusTitles),
             '' if len(self.virusTitles) == 1 else 'es',
             len(self.sampleNames),
             '' if len(self.sampleNames) == 1 else 's'))

    def maxProteinFraction(self, virusTitle):
        """
        Get the fraction of a virus' proteins matched by any sample that
        matches the virus.

        @param virusTitle: A C{str} virus title.
        @return: The C{float} maximum fraction of a virus' proteins that is
            matched by any sample. If the number of proteins for the virus is
            unknown, return 1.0 (i.e., assume all proteins are matched).
        """

        proteinCount = self._virusProteinCount[virusTitle]
        if proteinCount:
            maxMatches = max(
                len(sample['proteins'])
                for sample in self.virusTitles[virusTitle].values())
            return maxMatches / proteinCount
        else:
            return 1.0

    def addFile(self, filename, fp):
        """
        Read and record protein information for a sample.

        @param filename: A C{str} file name.
        @param fp: An open file pointer to read the file's data from.
        @raise ValueError: If information for a virus/protein/sample
            combination is given more than once.
        """

        if self._sampleNameRegex:
            match = self._sampleNameRegex.search(filename)
            if match:
                sampleName = match.group(1)
            else:
                sampleName = filename
        else:
            sampleName = filename

        outDir = join(dirname(filename), self._assetDir)

        self.sampleNames[sampleName] = join(outDir, 'index.html')

        for index, proteinLine in enumerate(fp):
            proteinLine = proteinLine[:-1]
            (coverage, medianScore, bestScore, readCount, hspCount,
             proteinLength, titles) = proteinLine.split(None, 6)

            proteinTitle, virusTitle = splitTitles(titles)

            if virusTitle not in self.virusTitles:
                self.virusTitles[virusTitle] = {}

            if sampleName not in self.virusTitles[virusTitle]:
                self.virusTitles[virusTitle][sampleName] = {
                    'proteins': {},
                    'uniqueReadCount': None,
                }

            proteins = self.virusTitles[virusTitle][sampleName]['proteins']

            # We should only receive one line of information for a given
            # virus/sample/protein combination.
            if proteinTitle in proteins:
                raise ValueError(
                    'Protein %r already seen for virus %r sample %r.' %
                    (proteinTitle, virusTitle, sampleName))

            proteins[proteinTitle] = {
                'bestScore': float(bestScore),
                'bluePlotFilename': join(outDir, '%d.png' % index),
                'coverage': float(coverage),
                'readsFilename': join(outDir, '%d.%s' % (index, self._format)),
                'hspCount': int(hspCount),
                'index': index,
                'medianScore': float(medianScore),
                'outDir': outDir,
                'proteinLength': int(proteinLength),
                'proteinTitle': proteinTitle,
                'proteinURL': NCBISequenceLinkURL(proteinTitle),
                'readCount': int(readCount),
            }

    def _computeUniqueReadCounts(self):
        """
        Add all virus / sample combinations to self.virusSampleFiles.

        This will make all de-duplicated FASTA/FASTQ files and store the number
        of de-duplicated reads into C{self.virusTitles}.
        """
        for virusTitle, samples in self.virusTitles.items():
            for sampleName in samples:
                self.virusSampleFiles.add(virusTitle, sampleName)

    def toStr(self):
        """
        Produce a string representation of the virus summary.

        @return: A C{str} suitable for printing.
        """
        # Note that the string representation contains much less
        # information than the HTML summary. E.g., it does not contain the
        # unique (de-duplicated) read count, since that is only computed
        # when we are making combined FASTA files of reads matching a
        # virus.
        readCountGetter = itemgetter('readCount')
        result = []
        append = result.append

        append(self._title())
        append('')

        for virusTitle in sorted(self.virusTitles):
            samples = self.virusTitles[virusTitle]
            sampleCount = len(samples)
            append('%s (in %d sample%s)' %
                   (virusTitle,
                    sampleCount, '' if sampleCount == 1 else 's'))
            for sampleName in sorted(samples):
                proteins = samples[sampleName]['proteins']
                proteinCount = len(proteins)
                totalReads = sum(readCountGetter(p) for p in proteins.values())
                append('  %s (%d protein%s, %d read%s)' %
                       (sampleName,
                        proteinCount, '' if proteinCount == 1 else 's',
                        totalReads, '' if totalReads == 1 else 's'))
                for proteinTitle in sorted(proteins):
                    append(
                        '    %(coverage).2f\t%(medianScore).2f\t'
                        '%(bestScore).2f\t%(readCount)4d\t%(hspCount)4d\t'
                        '%(index)3d\t%(proteinTitle)s'
                        % proteins[proteinTitle])
            append('')

        return '\n'.join(result)

    def toHTML(self, virusPanelFilename=None, minProteinFraction=0.0):
        """
        Produce an HTML string representation of the virus summary.

        @param virusPanelFilename: If not C{None}, a C{str} filename to write
            a virus panel PNG image to.
        @param minProteinFraction: The C{float} minimum fraction of proteins
            in a virus that must be matched by at least one sample in order for
            that virus to be displayed.
        @return: An HTML C{str} suitable for printing.
        """
        highlightSymbol = '&starf;'
        self._computeUniqueReadCounts()

        if virusPanelFilename:
            self.virusPanel(virusPanelFilename)

        virusTitles = sorted(
            virusTitle for virusTitle in self.virusTitles
            if self.maxProteinFraction(virusTitle) >= minProteinFraction)
        nVirusTitles = len(virusTitles)
        sampleNames = sorted(self.sampleNames)

        result = [
            '<html>',
            '<head>',
            '<title>',
            'Summary of viruses',
            '</title>',
            '<meta charset="UTF-8">',
            '</head>',
            '<body>',
            '<style>',
            '''\
            body {
                margin-left: 2%;
                margin-right: 2%;
            }
            hr {
                display: block;
                margin-top: 0.5em;
                margin-bottom: 0.5em;
                margin-left: auto;
                margin-right: auto;
                border-style: inset;
                border-width: 1px;
            }
            .significant {
                color: red;
                margin-right: 2px;
            }
            .sample {
                margin-bottom: 2px;
            }
            .indented {
                margin-left: 2em;
            }
            .sample-name {
                font-size: 125%;
                font-weight: bold;
            }
            .virus-title {
                font-size: 125%;
                font-weight: bold;
            }
            .index-title {
                font-weight: bold;
            }
            .index {
                font-size: small;
            }
            .protein-title {
                font-family: "Courier New", Courier, monospace;
            }
            .stats {
                font-family: "Courier New", Courier, monospace;
                white-space: pre;
            }
            .protein-list {
                margin-top: 2px;
            }''',
            '</style>',
            '</head>',
            '<body>',
        ]

        proteinFieldsDescription = (
            '<p>',
            'In all bullet point protein lists below, there are eight fields:',
            '<ol>',
            '<li>Coverage fraction.</li>',
            '<li>Median bit score.</li>',
            '<li>Best bit score.</li>',
            '<li>Read count.</li>',
            '<li>HSP count (a read can match a protein more than once).</li>',
            '<li>Protein length (in AAs).</li>',
            '<li>Index (just ignore this).</li>',
            '<li>Protein name.</li>',
            '</ol>',
            '</p>',
        )

        append = result.append

        append('<h1>Summary of viruses</h1>')
        append('<p>')
        append(self._title())

        if self._virusProteinCount:
            percent = minProteinFraction * 100.0
            if nVirusTitles < len(self.virusTitles):
                if nVirusTitles == 1:
                    append('Virus protein fraction filtering has been '
                           'applied, so information on only 1 virus is '
                           'displayed. This is the only virus for which at '
                           'least one sample matches at least %.2f%% of the '
                           'viral proteins.' % percent)
                else:
                    append('Virus protein fraction filtering has been '
                           'applied, so information on only %d viruses is '
                           'displayed. These are the only viruses for which '
                           'at least one sample matches at least %.2f%% of '
                           'the viral proteins.' % (nVirusTitles, percent))
            else:
                append('Virus protein fraction filtering has been applied, '
                       'but all viruses have at least %.2f%% of their '
                       'proteins matched by at least one sample.')

            append('Samples that match a virus (and viruses with a matching '
                   'sample) with at least this protein fraction are '
                   'highlighted using <span class="significant">%s</span>.' %
                   highlightSymbol)

        append('</p>')

        if virusPanelFilename:
            append('<p>')
            append('<a href="%s">Panel showing read count per virus, per '
                   'sample.</a>' % virusPanelFilename)
            append('Red vertical bars indicate samples with an unusually high '
                   'read count.')
            append('</p>')

        result.extend(proteinFieldsDescription)

        # Write a linked table of contents by virus.
        append('<p><span class="index-title">Virus index:</span>')
        append('<span class="index">')
        for virusTitle in virusTitles:
            append('<a href="#virus-%s">%s</a>' % (virusTitle, virusTitle))
            append('&middot;')
        # Get rid of final middle dot and add a period.
        result.pop()
        result[-1] += '.'
        append('</span></p>')

        # Write a linked table of contents by sample.
        append('<p><span class="index-title">Sample index:</span>')
        append('<span class="index">')
        for sampleName in sampleNames:
            append('<a href="#sample-%s">%s</a>' % (sampleName, sampleName))
            append('&middot;')
        # Get rid of final middle dot and add a period.
        result.pop()
        result[-1] += '.'
        append('</span></p>')

        # Write all viruses (with samples (with proteins)).
        append('<hr>')
        append('<h1>Viruses by sample</h1>')

        for virusTitle in virusTitles:
            samples = self.virusTitles[virusTitle]
            sampleCount = len(samples)
            virusProteinCount = self._virusProteinCount[virusTitle]
            append(
                '<a id="virus-%s"></a>'
                '<span class="virus-title"><a href="%s%s">%s</a></span>'
                '%s, was matched by %d sample%s:' %
                (virusTitle, self.VIRALZONE, quote(virusTitle),
                 virusTitle,
                 ((' (with %d protein%s)' %
                   (virusProteinCount, '' if virusProteinCount == 1 else 's'))
                  if virusProteinCount else ''),
                 sampleCount,
                 '' if sampleCount == 1 else 's'))
            for sampleName in sorted(samples):
                readsFileName = self.virusSampleFiles.lookup(virusTitle,
                                                             sampleName)
                proteins = samples[sampleName]['proteins']
                proteinCount = len(proteins)
                uniqueReadCount = samples[sampleName]['uniqueReadCount']
                if virusProteinCount and (
                        proteinCount / virusProteinCount >=
                        minProteinFraction):
                    highlight = ('<span class="significant">%s</span>' %
                                 highlightSymbol)
                else:
                    highlight = ''

                append(
                    '<p class="sample indented">'
                    '%sSample <a href="#sample-%s">%s</a> '
                    '(%d protein%s, <a href="%s">%d read%s</a>, '
                    '<a href="%s">panel</a>):</p>' %
                    (highlight, sampleName, sampleName,
                     proteinCount, '' if proteinCount == 1 else 's',
                     readsFileName,
                     uniqueReadCount, '' if uniqueReadCount == 1 else 's',
                     self.sampleNames[sampleName]))
                append('<ul class="protein-list indented">')
                for proteinTitle in sorted(proteins):
                    proteinMatch = proteins[proteinTitle]
                    append(
                        '<li>'
                        '<span class="stats">'
                        '%(coverage).2f %(medianScore).2f %(bestScore).2f '
                        '%(readCount)4d %(hspCount)4d %(proteinLength)4d '
                        '%(index)3d '
                        '</span> '
                        '<span class="protein-title">'
                        '%(proteinTitle)s'
                        '</span> '
                        '(<a href="%(bluePlotFilename)s">blue plot</a>, '
                        '<a href="%(readsFilename)s">reads</a>'
                        % proteinMatch)

                    if proteinMatch['proteinURL']:
                        # Append this directly to the last string in result, to
                        # avoid introducing whitespace when we join result
                        # using '\n'.
                        result[-1] += (', <a href="%s">NCBI</a>' %
                                       proteinMatch['proteinURL'])
                    result[-1] += ')'

                    append('</li>')

                append('</ul>')

        # Write all samples (with viruses (with proteins)).
        append('<hr>')
        append('<h1>Samples by virus</h1>')

        for sampleName in sampleNames:
            sampleVirusTitles = set()
            for virusTitle in virusTitles:
                if (sampleName in self.virusTitles[virusTitle] and
                        self.maxProteinFraction(virusTitle) >=
                        minProteinFraction):
                    sampleVirusTitles.add(virusTitle)

            append(
                '<a id="sample-%s"></a>'
                'Sample <span class="sample-name">%s</span> '
                'matched proteins from %d virus%s, '
                '<a href="%s">panel</a>:<br/>' %
                (sampleName, sampleName, len(sampleVirusTitles),
                 '' if len(sampleVirusTitles) == 1 else 'es',
                 self.sampleNames[sampleName]))

            for virusTitle in sorted(sampleVirusTitles):
                readsFileName = self.virusSampleFiles.lookup(virusTitle,
                                                             sampleName)
                proteins = self.virusTitles[virusTitle][sampleName]['proteins']
                uniqueReadCount = self.virusTitles[
                    virusTitle][sampleName]['uniqueReadCount']
                proteinCount = len(proteins)
                virusProteinCount = self._virusProteinCount[virusTitle]

                highlight = ''
                if virusProteinCount:
                    proteinCountStr = '%d/%d protein%s' % (
                        proteinCount, virusProteinCount,
                        '' if virusProteinCount == 1 else 's')
                    if proteinCount / virusProteinCount >= minProteinFraction:
                        highlight = ('<span class="significant">%s</span>' %
                                     highlightSymbol)
                else:
                    proteinCountStr = '%d protein%s' % (
                        proteinCount, '' if proteinCount == 1 else 's')

                append(
                    '<p class="sample indented">'
                    '%s<a href="#virus-%s">%s</a> %s, '
                    '<a href="%s">%d read%s</a>:</p>' %
                    (highlight, virusTitle, virusTitle,
                     proteinCountStr, readsFileName,
                     uniqueReadCount, '' if uniqueReadCount == 1 else 's'))
                append('<ul class="protein-list indented">')
                for proteinTitle in sorted(proteins):
                    proteinMatch = proteins[proteinTitle]
                    append(
                        '<li>'
                        '<span class="stats">'
                        '%(coverage).2f %(medianScore).2f %(bestScore).2f '
                        '%(readCount)4d %(hspCount)4d %(proteinLength)4d '
                        '%(index)3d '
                        '</span> '
                        '<span class="protein-title">'
                        '%(proteinTitle)s'
                        '</span> '
                        '(<a href="%(bluePlotFilename)s">blue plot</a>, '
                        '<a href="%(readsFilename)s">reads</a>'
                        % proteinMatch)

                    if proteinMatch['proteinURL']:
                        # Append this directly to the last string in result, to
                        # avoid introducing whitespace when we join result
                        # using '\n'.
                        result[-1] += (', <a href="%s">NCBI</a>' %
                                       proteinMatch['proteinURL'])
                    result[-1] += ')'

                    append('</li>')

                append('</ul>')

        append('</body>')
        append('</html>')

        return '\n'.join(result)

    def _virusSamplePlot(self, virusTitle, sampleNames, ax):
        """
        Make an image of a graph giving virus read count (Y axis) versus
        sample id (X axis).

        @param virusTitle: A C{str} virus title.
        @param sampleNames: A sorted C{list} of sample names.
        @param ax: A matplotlib C{axes} instance.
        """
        readCounts = []
        for i, sampleName in enumerate(sampleNames):
            try:
                readCount = self.virusTitles[virusTitle][sampleName][
                    'uniqueReadCount']
            except KeyError:
                readCount = 0
            readCounts.append(readCount)

        highlight = 'r'
        normal = 'gray'
        sdMultiple = 2.5
        minReadsForHighlighting = 10
        highlighted = []

        if len(readCounts) == 1:
            if readCounts[0] > minReadsForHighlighting:
                color = [highlight]
                highlighted.append(sampleNames[0])
            else:
                color = [normal]
        else:
            mean = np.mean(readCounts)
            sd = np.std(readCounts)
            color = []
            for readCount, sampleName in zip(readCounts, sampleNames):
                if (readCount > (sdMultiple * sd) + mean and
                        readCount >= minReadsForHighlighting):
                    color.append(highlight)
                    highlighted.append(sampleName)
                else:
                    color.append(normal)

        nSamples = len(sampleNames)
        x = np.arange(nSamples)
        yMin = np.zeros(nSamples)
        ax.set_xticks([])
        ax.set_xlim((-0.5, nSamples - 0.5))
        ax.vlines(x, yMin, readCounts, color=color)
        if highlighted:
            title = '%s\nIn red: %s' % (
                virusTitle, fill(', '.join(highlighted), 50))
        else:
            # Add a newline to keep the first line of each title at the
            # same place as those titles that have an "In red:" second
            # line.
            title = virusTitle + '\n'

        ax.set_title(title, fontsize=10)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.tick_params(axis='both', which='minor', labelsize=6)

    def virusPanel(self, filename):
        """
        Make a panel of images, with each image being a graph giving virus
        de-duplicated read count (Y axis) versus sample id (X axis).

        @param filename: A C{str} file name to write the image to.
        """
        self._computeUniqueReadCounts()
        virusTitles = sorted(self.virusTitles)
        sampleNames = sorted(self.sampleNames)

        cols = 5
        rows = int(len(virusTitles) / cols) + (
            0 if len(virusTitles) % cols == 0 else 1)
        figure, ax = plt.subplots(rows, cols, squeeze=False)

        coords = dimensionalIterator((rows, cols))

        for i, virusTitle in enumerate(virusTitles):
            row, col = next(coords)
            self._virusSamplePlot(virusTitle, sampleNames, ax[row][col])

        # Hide the final panel graphs (if any) that have no content. We do
        # this because the panel is a rectangular grid and some of the
        # plots at the end of the last row may be unused.
        for row, col in coords:
            ax[row][col].axis('off')

        figure.suptitle(
            ('Per-sample read count for %d virus%s and %d sample%s.\n\nSample '
             'name%s: %s') % (
                 len(virusTitles),
                 '' if len(virusTitles) == 1 else 'es',
                 len(sampleNames),
                 '' if len(sampleNames) == 1 else 's',
                 '' if len(sampleNames) == 1 else 's',
                 fill(', '.join(sampleNames), 50)),
            fontsize=20)
        figure.set_size_inches(5.0 * cols, 2.0 * rows, forward=True)
        plt.subplots_adjust(hspace=0.4)

        figure.savefig(filename)

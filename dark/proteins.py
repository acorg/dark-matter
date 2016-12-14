from __future__ import division, print_function

import re
from os.path import dirname, join
from operator import itemgetter
from six.moves.urllib.parse import quote
import numpy as np
from textwrap import fill

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
from dark.html import NCBISequenceLinkURL
from dark.reads import Reads


class VirusSampleFASTA(object):
    """
    Maintain a cache of virus/sample FASTA file names, creating de-duplicaed
    FASTA files (from reads for all proteins of a virus that a sample has)
    on demand.

    @param proteinGrouper: An instance of C{ProteinGrouper}.
    """
    def __init__(self, proteinGrouper):
        self._proteinGrouper = proteinGrouper
        self._viruses = {}
        self._samples = {}
        self._fastaFilenames = {}

    def add(self, virusTitle, sampleName):
        """
        Add a virus title, sample name combination and get its FASTA file
        name and unique read count. Write the FASTA file if it does not
        already exist. Save the unique read count into
        C{self._proteinGrouper}.

        @param virusTitle: A C{str} virus title.
        @param sampleName: A C{str} sample name.
        @return: A C{str} giving the FASTA file name holding all the reads
            (without duplicates) from the sample that matched the proteins in
            the given virus.
        """
        virusIndex = self._viruses.setdefault(virusTitle, len(self._viruses))
        sampleIndex = self._samples.setdefault(sampleName, len(self._samples))

        try:
            return self._fastaFilenames[(virusIndex, sampleIndex)]
        except KeyError:
            reads = Reads()
            for protein in self._proteinGrouper.virusTitles[
                    virusTitle][sampleName]['proteins']:
                for read in FastaReads(protein['fastaFilename'],
                                       checkAlphabet=0):
                    reads.add(read)
            saveFilename = join(
                protein['outDir'],
                'virus-%d-sample-%d.fasta' % (virusIndex, sampleIndex))
            reads.filter(removeDuplicates=True)
            nReads = reads.save(saveFilename)
            # Save the unique read count into self._proteinGrouper
            self._proteinGrouper.virusTitles[
                virusTitle][sampleName]['uniqueReadCount'] = nReads
            self._fastaFilenames[(virusIndex, sampleIndex)] = saveFilename
            return saveFilename

    def lookup(self, virusTitle, sampleName):
        """
        Look up a virus title, sample name combination and get its FASTA file
        name and unique read count.

        This method should be used instead of C{add} in situations where
        you want an exception to be raised if a virus/sample combination has
        not already been passed to C{add}.

        @param virusTitle: A C{str} virus title.
        @param sampleName: A C{str} sample name.
        @raise KeyError: If the virus title or sample name have not been seen,
            either individually or in combination.
        @return: A (C{str}, C{int}) tuple retrieved from self._fastaFilenames
        """
        virusIndex = self._viruses[virusTitle]
        sampleIndex = self._samples[sampleName]
        return self._fastaFilenames[(virusIndex, sampleIndex)]


class ProteinGrouper(object):
    """
    Group matched proteins by the virus they come from.

    @param assetDir: The C{str} directory name where
        C{noninteractive-alignment-panel.py} put its HTML, blue plot and
        alignment panel images, and FASTA files. This must be relative
        to the filenames that will later be passed to C{addFile}.
    @param sampleNameRegex: A C{str} regular expression that can be used to
        extract a short sample name from full file names subsequently passed
        to C{self.addFile}. The regular expression must have a matching group
        (delimited by parentheses) to capture the part of the file name that
        should be used as the sample name.
    """

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
    # Unfortunately the regex doesn't find the vius title when the protein
    # title has nested [...] sections, as in this example:
    #
    #   gi|224808893|ref|YP_002643049.1| replication-associated protein [Tomato
    #   leaf curl Nigeria virus-[Nigeria:2006]]
    #
    # I decided not to worry about nested [...] sections (there are only 2
    # instances that I know of).
    VIRUS_RE = re.compile('^(.*)\[([^\]]+)\]$')

    # The virus title assigned to proteins whose title strings cannot be parsed
    # for a virus title (see previous comment).  Do not use '<', '>' or any
    # other HTML special chars in the following.
    NO_VIRUS_TITLE = '[no virus name found in protein title]'

    VIRALZONE = 'http://viralzone.expasy.org/cgi-bin/viralzone/search?query='

    def __init__(self, assetDir='out', sampleNameRegex=None):
        self._sampleNameRegex = (re.compile(sampleNameRegex) if sampleNameRegex
                                 else None)
        self._assetDir = assetDir
        # virusTitles will be a dict of dicts of dicts. The first two keys
        # will be a virus title and a sample name. The final dict will
        # contain 'proteins' (a list of dicts) and 'uniqueReadCount' (an int).
        self.virusTitles = {}
        # sampleNames is keyed by sample name and will have values that hold
        # the sample's alignment panel index.html file.
        self.sampleNames = {}
        self.virusSampleFASTA = VirusSampleFASTA(self)

    def _title(self):
        """
        Create a title summarizing the viruses and samples.

        @return: A C{str} title.
        """
        return (
            '%d virus%s found in %d sample%s' %
            (len(self.virusTitles), '' if len(self.virusTitles) == 1 else 'es',
             len(self.sampleNames), '' if len(self.sampleNames) == 1 else 's'))

    def addFile(self, filename, fp):
        """
        Read and record protein information for a sample.

        @param filename: A C{str} file name.
        @param fp: An open file pointer to read the file's data from.
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

            match = self.VIRUS_RE.match(titles)
            if match:
                proteinTitle = match.group(1).strip()
                virusTitle = match.group(2)
            else:
                proteinTitle = titles
                virusTitle = self.NO_VIRUS_TITLE

            if virusTitle not in self.virusTitles:
                self.virusTitles[virusTitle] = {}

            if sampleName not in self.virusTitles[virusTitle]:
                self.virusTitles[virusTitle][sampleName] = {
                    'proteins': [],
                    'uniqueReadCount': None,
                }

            self.virusTitles[virusTitle][sampleName]['proteins'].append({
                'bestScore': float(bestScore),
                'bluePlotFilename': join(outDir, '%d.png' % index),
                'coverage': float(coverage),
                'fastaFilename': join(outDir, '%d.fasta' % index),
                'hspCount': int(hspCount),
                'index': index,
                'medianScore': float(medianScore),
                'outDir': outDir,
                'proteinLength': int(proteinLength),
                'proteinTitle': proteinTitle,
                'proteinURL': NCBISequenceLinkURL(proteinTitle),
                'readCount': int(readCount),
            })

    def _computeUniqueReadCounts(self):
        """
        Add all virus / sample combinations to self.virusSampleFASTA.

        This will make all de-duplicated FASTA files and store the number
        of de-duplicated reads into C{self.virusTitles}.
        """
        for virusTitle, samples in self.virusTitles.items():
            for sampleName in samples:
                self.virusSampleFASTA.add(virusTitle, sampleName)

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
        titleGetter = itemgetter('proteinTitle')
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
                totalReads = sum(readCountGetter(p) for p in proteins)
                append('  %s (%d protein%s, %d read%s)' %
                       (sampleName,
                        proteinCount, '' if proteinCount == 1 else 's',
                        totalReads, '' if totalReads == 1 else 's'))
                proteins.sort(key=titleGetter)
                for proteinMatch in proteins:
                    append(
                        '    %(coverage).2f\t%(medianScore).2f\t'
                        '%(bestScore).2f\t%(readCount)4d\t%(hspCount)4d\t'
                        '%(index)3d\t%(proteinTitle)s'
                        % proteinMatch)
            append('')

        return '\n'.join(result)

    def toHTML(self, virusPanelFilename=None):
        """
        Produce an HTML string representation of the virus summary.

        @param virusPanelFilename: If not C{None}, a C{str} filename to write
            a virus panel PNG image to.
        @return: An HTML C{str} suitable for printing.
        """
        self._computeUniqueReadCounts()

        if virusPanelFilename:
            self.virusPanel(virusPanelFilename)

        titleGetter = itemgetter('proteinTitle')
        virusTitles = sorted(self.virusTitles)
        sampleNames = sorted(self.sampleNames)

        result = [
            '<html>',
            '<head>',
            '<title>',
            self._title(),
            '</title>',
            '</head>',
            '<body>',
            '<style>',
            '''\
            body {
                margin-left: 2%;
                margin-right: 2%;
            }
            .sample {
                margin-bottom: 2px;
            }
            .sample-name {
                color: red;
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
            'In the bullet point protein lists below, there are eight fields:',
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

        append('<h1>%s</h1>' % self._title())

        if virusPanelFilename:
            append('<p>')
            append('<a href="%s">Panel showing read count per virus, per '
                   'sample.</a>' % virusPanelFilename)
            append('Red vertical bars indicate samples with an unusually high '
                   'read count.')
            append('</p>')

        # Write a linked table of contents by virus.
        append('<h2>Virus index</h2>')
        append('<p class="index">')
        for virusTitle in virusTitles:
            append('<a href="#virus-%s">%s</a>' % (virusTitle, virusTitle))
            append('&middot;')
        # Get rid of final middle dot.
        result.pop()
        append('</p>')

        # Write a linked table of contents by sample.
        append('<h2>Sample index</h2>')
        append('<p class="index">')
        for sampleName in sampleNames:
            append('<a href="#sample-%s">%s</a>' % (sampleName, sampleName))
            append('&middot;')
        # Get rid of final middle dot.
        result.pop()
        append('</p>')

        # Write all viruses (with samples (with proteins)).
        append('<h1>Viruses by sample</h1>')
        result.extend(proteinFieldsDescription)

        for virusTitle in virusTitles:
            samples = self.virusTitles[virusTitle]
            sampleCount = len(samples)
            append(
                '<a id="virus-%s">'
                '<h2 class="virus"><span class="virus-title">%s</span> '
                '(in %d sample%s, <a href="%s%s">viralzone</a>)</h2>' %
                (virusTitle, virusTitle, sampleCount,
                 '' if sampleCount == 1 else 's',
                 self.VIRALZONE, quote(virusTitle)))
            for sampleName in sorted(samples):
                fastaName = self.virusSampleFASTA.lookup(virusTitle,
                                                         sampleName)
                proteins = samples[sampleName]['proteins']
                proteinCount = len(proteins)
                uniqueReadCount = samples[sampleName]['uniqueReadCount']
                append(
                    '<p class=sample>'
                    '<span class="sample-name">%s</span> '
                    '(%d protein%s, %d read%s, <a href="%s">panel</a>, '
                    '<a href="%s">fasta</a>)' %
                    (sampleName,
                     proteinCount, '' if proteinCount == 1 else 's',
                     uniqueReadCount, '' if uniqueReadCount == 1 else 's',
                     self.sampleNames[sampleName], fastaName))
                proteins.sort(key=titleGetter)
                append('<ul class="protein-list">')
                for proteinMatch in proteins:
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
                        '<a href="%(fastaFilename)s">fasta</a>'
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
                append('</p>')

        # Write all samples (with viruses (with proteins)).
        append('<h1>Samples by virus</h1>')
        result.extend(proteinFieldsDescription)

        for sampleName in sampleNames:
            sampleVirusTitles = set()
            for virusTitle in virusTitles:
                if sampleName in self.virusTitles[virusTitle]:
                    sampleVirusTitles.add(virusTitle)

            append(
                '<a id="sample-%s">'
                '<h2 class="sample"><span class="sample-name">%s</span> '
                '(has proteins from %d virus%s, <a href="%s">panel</a>)</h2>' %
                (sampleName, sampleName, len(sampleVirusTitles),
                 '' if len(sampleVirusTitles) == 1 else 'es',
                 self.sampleNames[sampleName]))

            for virusTitle in sorted(sampleVirusTitles):
                fastaName = self.virusSampleFASTA.lookup(virusTitle,
                                                         sampleName)
                proteins = self.virusTitles[virusTitle][sampleName]['proteins']
                uniqueReadCount = self.virusTitles[
                    virusTitle][sampleName]['uniqueReadCount']
                proteinCount = len(proteins)
                append(
                    '<p class="sample">'
                    '<span class="virus-title">%s</span> '
                    '(%d protein%s, %d read%s, <a href="%s">fasta</a>)' %
                    (virusTitle,
                     proteinCount, '' if proteinCount == 1 else 's',
                     uniqueReadCount, '' if uniqueReadCount == 1 else 's',
                     fastaName))
                proteins.sort(key=titleGetter)
                append('<ul class="protein-list">')
                for proteinMatch in proteins:
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
                        '<a href="%(fastaFilename)s">fasta</a>'
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
                append('</p>')

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

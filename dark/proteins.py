import re
from os.path import dirname, join
from collections import defaultdict
from operator import itemgetter
from six.moves.urllib.parse import quote

from dark.html import NCBISequenceLinkURL


class ProteinGrouper(object):
    """
    Group matched proteins by the virus they come from.

    @param sampleNameRegex: A C{str} that ......
    """

    # The following regex is deliberately greedy (using .*) to consume the
    # whole title before backtracking. That way it will match just the last
    # [title] in a protein. This avoids situations in which two [...] delimited
    # substrings are present in a protein name (in which case we just want the
    # last). E.g., the following is a complete protein name:
    #
    #   gi|19919894|ref|NP_612577.1| Enzymatic polyprotein [Contains: Aspartic
    #   protease; Endonuclease; Reverse transcriptase] [Carnation etched ring
    #   virus]
    #
    # Unfortunately this doesn't match names like
    #
    #   gi|224808893|ref|YP_002643049.1| replication-associated protein [Tomato
    #   leaf curl Nigeria virus-[Nigeria:2006]]
    #
    # but I decided not to worry about them (there are only 2 instances that I
    # know of).
    VIRUS_RE = re.compile('^(.*)\[([^\]]+)\]$')

    # Do not use '<', '>' or any other HTML special chars in the following.
    NO_VIRUS_TITLE = '[no virus name found in protein title]'

    VIRALZONE = 'http://viralzone.expasy.org/cgi-bin/viralzone/search?query='

    def __init__(self, sampleNameRegex=None):
        self._sampleNameRegex = (re.compile(sampleNameRegex) if sampleNameRegex
                                 else None)
        # virusTitles is a dict of dicts of lists.
        self.virusTitles = defaultdict(lambda: defaultdict(list))
        # sampleNames has values that give the sample panel's index.html file.
        self.sampleNames = {}

    def addFile(self, filename, fp):
        """
        Read protein input and add it to virusTitles.

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

        # outDir is the assumed location of the directory where
        # noninteractive-alignment-panel.py puts its HTML, image, and FASTA
        # output.
        outDir = join(dirname(filename), 'out')

        self.sampleNames[sampleName] = join(outDir, 'index.html')

        for index, proteinLine in enumerate(fp):
            proteinLine = proteinLine[:-1]
            (coverage, medianScore, bestScore, readCount, hspCount,
             virusLength, titles) = proteinLine.split('\t')

            match = self.VIRUS_RE.match(titles)
            if match:
                proteinTitle = match.group(1).strip()
                virusTitle = match.group(2)
            else:
                proteinTitle = titles
                virusTitle = self.NO_VIRUS_TITLE

            self.virusTitles[virusTitle][sampleName].append({
                'bluePlotFilename': join(outDir, '%d.png' % index),
                'coverage': float(coverage),
                'fastaFilename': join(outDir, '%d.fasta' % index),
                'index': index,
                'medianScore': float(medianScore),
                'bestScore': float(bestScore),
                'readCount': int(readCount),
                'hspCount': int(hspCount),
                'proteinTitle': proteinTitle,
                'proteinURL': NCBISequenceLinkURL(proteinTitle),
            })

    def toStr(self):
        """
        """
        titleGetter = itemgetter('proteinTitle')
        readCountGetter = itemgetter('readCount')
        result = []
        append = result.append

        append('%d viruses found in %d samples' %
               (len(self.virusTitles), len(self.sampleNames)))
        append('')

        for virusTitle in sorted(self.virusTitles):
            samples = self.virusTitles[virusTitle]
            sampleCount = len(samples)
            append('%s (in %d sample%s)' %
                   (virusTitle,
                    sampleCount, '' if sampleCount == 1 else 's'))
            for sampleName in sorted(samples):
                proteins = samples[sampleName]
                proteinCount = len(proteins)
                totalReads = sum(readCountGetter(p) for p in proteins)
                append('  %s (%d protein%s, %d read%s)' %
                       (sampleName,
                        proteinCount, '' if proteinCount == 1 else 's',
                        totalReads, '' if totalReads == 1 else 's'))
                proteins.sort(key=titleGetter)
                for proteinMatch in proteins:
                    append(
                        '    %(coverage)f\t%(medianScore)f\t%(bestScore)f\t'
                        '%(readCount)d\t%(hspCount)d\t%(index)d\t'
                        '%(proteinTitle)s'
                        % proteinMatch)
            append('')

        return '\n'.join(result)

    def toHTML(self):
        """
        """
        titleGetter = itemgetter('proteinTitle')
        readCountGetter = itemgetter('readCount')
        result = [
            '<html>',
            '<head>',
            '<title>',
            '%d viruses found in %d samples' % (len(self.virusTitles),
                                                len(self.sampleNames)),
            '</title>',
            '</head>',
            '<body>',
            '<style>',
            '''
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
            }
            ''',
            '</style>',
            '</head>',
            '<body>',
        ]

        append = result.append

        virusTitles = sorted(self.virusTitles)
        sampleNames = sorted(self.sampleNames)

        append('<h1>%d viruses found in %d samples</h1>' %
               (len(virusTitles), len(sampleNames)))

        # Write a table of contents by virus.
        append('<h2>Virus index</h2>')
        append('<p class="index">')
        for index, virusTitle in enumerate(virusTitles):
            append('<a href="#virus-%d">%s</a>' % (index, virusTitle))
            append('&middot;')
        # Get rid of final middle dot.
        result.pop()
        append('</p>')

        # Write a table of contents by sample.
        append('<h2>Sample index</h2>')
        append('<p class="index">')
        for index, sampleName in enumerate(sampleNames):
            append('<a href="#sample-%d">%s</a>' % (index, sampleName))
            append('&middot;')
        # Get rid of final middle dot.
        result.pop()
        append('</p>')

        # Write all viruses (with samples (with proteins)).
        append('<h1>Viruses by sample</h1>')
        for index, virusTitle in enumerate(virusTitles):
            samples = self.virusTitles[virusTitle]
            sampleCount = len(samples)
            append(
                '<a id="virus-%d">'
                '<h2 class="virus"><span class="virus-title">%s</span> '
                '(in %d sample%s, <a href="%s%s">viralzone</a>)</h2>' %
                (index, virusTitle, sampleCount,
                 '' if sampleCount == 1 else 's',
                 self.VIRALZONE, quote(virusTitle)))
            for sampleName in sorted(samples):
                proteins = samples[sampleName]
                proteinCount = len(proteins)
                totalReads = sum(readCountGetter(p) for p in proteins)
                append(
                    '<p class=sample>'
                    '<span class="sample-name">%s</span> '
                    '(%d protein%s, %d read%s, <a href="%s">panel</a>)' %
                    (sampleName,
                     proteinCount, '' if proteinCount == 1 else 's',
                     totalReads, '' if totalReads == 1 else 's',
                     self.sampleNames[sampleName]))
                proteins.sort(key=titleGetter)
                append('<ul class="protein-list">')
                for proteinMatch in proteins:
                    append(
                        '<li>'
                        '<span class="stats">'
                        '%(coverage).2f %(medianScore).2f %(bestScore).2f '
                        '%(readCount)3d %(hspCount)3d %(index)3d '
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
        for index, sampleName in enumerate(sorted(sampleNames)):

            sampleVirusTitles = set()
            for virusTitle in virusTitles:
                if sampleName in self.virusTitles[virusTitle]:
                    sampleVirusTitles.add(virusTitle)

            append(
                '<a id="sample-%d">'
                '<h2 class="sample"><span class="sample-name">%s</span> '
                '(has proteins from %d virus%s, <a href="%s">panel</a>)</h2>' %
                (index, sampleName, len(sampleVirusTitles),
                 '' if len(sampleVirusTitles) == 1 else 'es',
                 self.sampleNames[sampleName]))

            for virusTitle in sorted(sampleVirusTitles):
                proteins = self.virusTitles[virusTitle][sampleName]
                proteinCount = len(proteins)
                totalReads = sum(readCountGetter(p) for p in proteins)
                append(
                    '<p class="sample">'
                    '<span class="virus-title">%s</span> '
                    '(%d protein%s, %d read%s)' %
                    (virusTitle,
                     proteinCount, '' if proteinCount == 1 else 's',
                     totalReads, '' if totalReads == 1 else 's'))
                proteins.sort(key=titleGetter)
                append('<ul class="protein-list">')
                for proteinMatch in proteins:
                    append(
                        '<li>'
                        '<span class="stats">'
                        '%(coverage).2f %(medianScore).2f %(bestScore).2f '
                        '%(readCount)3d %(hspCount)3d %(index)3d '
                        '</span> '
                        '<span class="protein-title">'
                        '%(proteinTitle)s'
                        '</span>'
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

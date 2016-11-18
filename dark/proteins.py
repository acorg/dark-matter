import re
from collections import defaultdict
from operator import itemgetter


class ProteinGrouper(object):
    """
    Group matched proteins by the virus they come from.

    @param sampleNameRegex: A C{str} that ......
    """

    VIRUS_RE = re.compile('^([^\[]+)\[([^\]]+)\]$')
    NO_VIRUS_TITLE = '<untitled>'

    def __init__(self, sampleNameRegex=None):
        self._sampleNameRegex = (re.compile(sampleNameRegex) if sampleNameRegex
                                 else None)
        # virusTitles is a dict of dicts of lists.
        self.virusTitles = defaultdict(lambda: defaultdict(list))

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
                'coverage': float(coverage),
                'index': index,
                'medianScore': float(medianScore),
                'bestScore': float(bestScore),
                'readCount': int(readCount),
                'hspCount': int(hspCount),
                'proteinTitle': proteinTitle,
            })

    def toStr(self):
        """
        """
        titleGetter = itemgetter('proteinTitle')
        readCountGetter = itemgetter('readCount')
        result = []
        append = result.append

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
            'Proteins grouped into viruses',
            '</title>',
            '</head>',
            '<body>',
            '<style>',
            '''
            .sample {
                /* color: green; */
            }
            .sample-name {
                color: red;
            }
            '''
            '</style>',
        ]

        append = result.append

        for virusTitle in sorted(self.virusTitles):
            samples = self.virusTitles[virusTitle]
            sampleCount = len(samples)
            append(
                '<h1 class="virus"><span class="virus-name">%s</span> '
                '(in %d sample%s)</h1>' %
                (virusTitle, sampleCount, '' if sampleCount == 1 else 's'))
            for sampleName in sorted(samples):
                proteins = samples[sampleName]
                proteinCount = len(proteins)
                totalReads = sum(readCountGetter(p) for p in proteins)
                append(
                    '<p class=sample>'
                    '<span class="sample-name">%s</span> '
                    '(%d protein%s, %d read%s)' %
                    (sampleName,
                     proteinCount, '' if proteinCount == 1 else 's',
                     totalReads, '' if totalReads == 1 else 's'))
                proteins.sort(key=titleGetter)
                append('<ul class="proteins">')
                for proteinMatch in proteins:
                    append(
                        '<li>'
                        '%(coverage)f\t%(medianScore)f\t%(bestScore)f\t'
                        '%(readCount)d\t%(hspCount)d\t%(index)d\t'
                        '%(proteinTitle)s'
                        '</li>'
                        % proteinMatch)
                append('</ul>')
                append('</p>')
            append('')
        append('</body>')
        append('</html>')

        return '\n'.join(result)

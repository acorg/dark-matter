from __future__ import print_function

from os.path import join

from dark.fastq import FastqReads
from dark.html import NCBISequenceLinkURL


class AlignmentPanelHTMLWriter(object):
    """
    Produces HTML details of a rectangular panel of graphs that each
    contain an alignment graph against a given sequence. This is
    supplementary output info for the AlignmentPanel class in graphics.py.

    @param outputDir: The C{str} directory to write files into.
    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    """
    def __init__(self, outputDir, titlesAlignments):
        self._outputDir = outputDir
        self._titlesAlignments = titlesAlignments
        self._images = []

    def addImage(self, imageBasename, accession, title, graphInfo):
        self._images.append({
            'accession': accession,
            'graphInfo': graphInfo,
            'imageBasename': imageBasename,
            'title': title
        })

    def close(self):
        with open(join(self._outputDir, 'index.html'), 'w') as fp:
            self._writeHeader(fp)
            self._writeBody(fp)
            self._writeFooter(fp)
        with open(join(self._outputDir, 'style.css'), 'w') as fp:
            self._writeCSS(fp)

    def _writeHeader(self, fp):
        fp.write("""\
<html>
  <head>
    <title>Read alignments for %d matched subjects</title>
    <link rel="stylesheet" type="text/css" href="style.css">
  </head>
  <body>
    <div id="content">
        """ % len(self._images))

    def _writeBody(self, fp):
        fp.write('<h1>Read alignments for %d matched subjects</h1>\n' %
                 len(self._images))

        # Write out an alignment panel as a table.
        cols = 6
        fp.write('<table><tbody>\n')

        for i, image in enumerate(self._images):
            title = image['title']
            accession = image['accession']
            if i % cols == 0:
                fp.write('<tr>\n')

            fp.write(
                '<td><a id="small_%s"></a><a href="#big_%s"><img src="%s" '
                'class="thumbnail"/></a></td>\n' %
                (accession, accession, image['imageBasename']))

            if i % cols == cols - 1:
                fp.write('</tr>')

        # Add empty cells to the final table row, and close the row, if
        # necessary.
        if i % cols < cols - 1:
            while i % cols < cols - 1:
                fp.write('<td>&nbsp;</td>\n')
                i += 1
            fp.write('</tr>\n')

        fp.write('</tbody></table>\n')

        # Write out the full images with additional detail.
        for i, image in enumerate(self._images):
            title = image['title']
            accession = image['accession']
            titleAlignments = self._titlesAlignments[title]
            graphInfo = image['graphInfo']
            readFormat = self._writeReads(image)
            fp.write("""
      <a id="big_%s"></a>
      <h3>%d: %s</h3>
      <p>
            Length: %d.
            Read count: %d.
            HSP count: %d.
            <a href="%s.%s">%s</a>.
            <a href="#small_%s">Top panel.</a>
"""
                     % (accession,
                        i, title,
                        titleAlignments.subjectLength,
                        titleAlignments.readCount(),
                        titleAlignments.hspCount(),
                        accession, readFormat, readFormat,
                        accession))

            url = NCBISequenceLinkURL(title)
            if url:
                fp.write('<a href="%s" target="_blank">NCBI</a>.' % url)

            # Write out feature information.
            if graphInfo['features'] is None:
                # Feature lookup was False (or we were offline).
                pass
            elif len(graphInfo['features']) == 0:
                fp.write('There were no features.')
            else:
                fp.write('<a href="%s">Features</a>' %
                         self._writeFeatures(i, image))

            # Write out the titles that this title invalidated due to its
            # read set.
            readSetFilter = self._titlesAlignments.readSetFilter
            if readSetFilter:
                invalidated = readSetFilter.invalidates(title)
                if invalidated:
                    nInvalidated = len(invalidated)
                    fp.write('<br/>This title invalidated %d other%s due to '
                             'its read set:<ul>'
                             % (nInvalidated,
                                '' if nInvalidated == 1 else 's'))
                    for title in invalidated:
                        fp.write('<li>%s</li>' % title)
                    fp.write('</ul>')

            fp.write('</p><img src="%s" class="full-size"/>' %
                     image['imageBasename'])

    def _writeFooter(self, fp):
        fp.write("""\
    </div>
  </body>
</html>
""")

    def _writeCSS(self, fp):
        fp.write("""\
#content {
  width: 95%;
  margin: auto;
}
img.thumbnail {
  height: 300px;
}
img.full-size {
  height: 900px;
}
""")

    def _writeReads(self, image):
        """
        Write a FASTA or FASTQ file containing the set of reads that hit a
        sequence.

        @param image: A member of self._images.
        @return: A C{str}, either 'fasta' or 'fastq' indicating the format
            of the reads in C{self._titlesAlignments}.
        """
        if isinstance(self._titlesAlignments.readsAlignments.reads,
                      FastqReads):
            format_ = 'fastq'
        else:
            format_ = 'fasta'
        filename = join(self._outputDir,
                        '%s.%s' % (image['accession'], format_))
        titleAlignments = self._titlesAlignments[image['title']]
        with open(filename, 'w') as fp:
            for titleAlignment in titleAlignments:
                fp.write(titleAlignment.read.toString(format_))
        return format_

    def _writeFeatures(self, image):
        """
        Write a text file containing the features as a table.

        @param image: A member of self._images.
        @return: The C{str} features file name - just the base name, not
            including the path to the file.
        """
        basename = image['accession'] + '-features.txt'
        filename = join(self._outputDir, basename)
        featureList = image['graphInfo']['features']
        # Note that the following (deliberately) creates an empty features
        # file if there were no features.
        with open(filename, 'w') as fp:
            for feature in featureList:
                fp.write('%s\n\n' % feature.feature)
        return basename

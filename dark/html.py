from Bio import SeqIO


def NCBISequenceLinkURL(title, default=None):
    """
    Given a sequence title, like "gi|42768646|gb|AY516849.1| Homo sapiens",
    return the URL of a link to the info page at NCBI.

    title: the sequence title to produce a link URL for.
    default: the value to return if the title cannot be parsed.
    """
    try:
        ref = title.split('|')[3].split('.')[0]
    except IndexError:
        return default
    else:
        return 'http://www.ncbi.nlm.nih.gov/nuccore/%s' % (ref,)


def NCBISequenceLink(title, default=None):
    """
    Given a sequence title, like "gi|42768646|gb|AY516849.1| Homo sapiens",
    return an HTML A tag dispalying a link to the info page at NCBI.

    title: the sequence title to produce an HTML link for.
    default: the value to return if the title cannot be parsed.
    """
    url = NCBISequenceLinkURL(title)
    if url is None:
        return default
    else:
        return '<a href="%s" target="_blank">%s</a>' % (url, title)


class AlignmentPanelHTML(object):
    """
    Produces HTML details of a rectangular panel of graphs that each
    contain an alignment graph against a given sequence. This is
    supplementary output info for the AlignmentPanel class in utils.py.

    @param outputDir: The C{str} directory to write files into.
    @param blastHits: A L{BlastHits} instance.
    """
    def __init__(self, outputDir, blastHits):
        self._outputDir = outputDir
        self._blastHits = blastHits
        self._images = []

    def addImage(self, imageBasename, title, alignmentInfo, plotInfo):
        self._images.append({
            'alignmentInfo': alignmentInfo,
            'plotInfo': plotInfo,
            'imageBasename': imageBasename,
            'title': title
        })

    def close(self, panelFilename):
        with open('%s/index.html' % self._outputDir, 'w') as fp:
            self._writeHeader(fp)
            self._writeBody(fp, panelFilename)
            self._writeFooter(fp)
        with open('%s/style.css' % self._outputDir, 'w') as fp:
            self._writeCSS(fp)

    def _writeHeader(self, fp):
        fp.write("""\
<html>
  <head>
    <title>woot</title>
    <link rel="stylesheet" type="text/css" href="style.css">
  </head>
  <body>
    <div id="content">
""")

    def _writeBody(self, fp, panelFilename):
        fp.write("""
      <h1>
        Original panel image
      </h1>
      <p>
        <img src="%s" class="panel"/>
      </p>
"""
                 % panelFilename)

        # Write out the summary images.
        for i, image in enumerate(self._images):
            title = image['title']
            alignmentInfo = image['alignmentInfo']
            plotInfo = image['plotInfo']
            readIds = self._writeFASTA(i, image)
            fp.write("""
      <a id="small_%d"></a>
      <h3>%d: %s</h3>
      <p>
        <a href="#big_%d"><img src="%s" class="thumbnail"/></a>
        Sequence length: %d base pairs.<br/>
        Number of reads that hit overall: %d.<br/>
        Number of HSPs in this selection: %d (of %d overall).<br/>
        <a href="#big_%d">Full size image</a>.
"""
                     % (i, i, title, i, image['imageBasename'],
                        self._blastHits.titles[title]['length'],
                        self._blastHits.titles[title]['readCount'],
                        len(plotInfo['items']), plotInfo['hspTotal'], i))

            url = NCBISequenceLinkURL(title)
            if url:
                fp.write('<br/><a href="%s" target="_blank">NCBI info on this '
                         'target</a>.' % url)

            # Write out feature information.
            if alignmentInfo['features'] is None:
                fp.write('<br/>Feature lookup was False (or we were offline).')
            elif len(alignmentInfo['features']) == 0:
                fp.write('<br/>There were no features.')
            else:
                fp.write('<br/><a href="%s">Features</a>' %
                         self._writeFeatures(i, image))

            # Write out the titles that this title invalidated due to its
            # read set.
            if self._blastHits.readSetFilter:
                invalidated = self._blastHits.readSetFilter.invalidates(title)
                if invalidated:
                    nInvalidated = len(invalidated)
                    fp.write('<br/>This title invalidated %d other%s due to '
                             'its read set:<ul>'
                             % (nInvalidated,
                                '' if nInvalidated == 1 else 's'))
                    for title in invalidated:
                        fp.write('<li>%s</li>' % title)
                    fp.write('</ul>')

            if len(plotInfo['items']):
                fp.write('<br/>Reads: <span class="reads">%s</span>'
                         % ', '.join(readIds))

            fp.write('<br clear="all"/></p>')

        # Write out the large images.
        for i, image in enumerate(self._images):
            title = image['title']
            fp.write("""
       <a id="big_%d"></a>
       <h3>%d: %s</h3>
       <p>
         <img src="%s"/>
         <a href="#small_%d">Back to summary</a>.<br/>
         <br clear="all"/>
       </p>
"""
                     % (i, i, title, image['imageBasename'], i))

    def _writeFooter(self, fp):
        fp.write("""\
    </div>
  </body>
</html>
""")

    def _writeCSS(self, fp):
        fp.write("""\
#content {
  width: 90%;
  margin: auto;
}
img.thumbnail {
  height: 200px;
  float: left;
}
img.panel {
  max-width: 1200px;
}
span.reads {
  font-size: 8pt;
  font-family: monospace;
}
""")

    def _writeFASTA(self, i, image):
        """
        Write a FASTA file containing the set of reads that hit a sequence.
        Return a list of the read identifiers.

        i: The number of the image in self._images.
        image: A member of self._images.
        """
        ids = []
        reads = []
        with open('%s/%d.fasta' % (self._outputDir, i), 'w') as fp:
            for item in image['plotInfo']['items']:
                readNum = item['readNum']
                reads.append(self._blastHits.fasta[readNum])
                ids.append(self._blastHits.fasta[readNum].description)
            SeqIO.write(reads, fp, 'fasta')
        return ids

    def _writeFeatures(self, i, image):
        """
        Write a txt file containing the features as a table.

        @param i: The number of the image in self._images.
        @param image: A member of self._images.
        @return: The C{str} features file name.
        """
        filename = '%s/features-%d.txt' % (self._outputDir, i)
        featureList = image['alignmentInfo']['features']
        with open(filename, 'w') as fp:
            for feature in featureList:
                fp.write(feature.description() + '\n')
        return filename

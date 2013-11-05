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

    fasta: A list of fasta sequences from SeqIO.parse
    outputDir: The directory to write files into.

    """

    def __init__(self, outputDir, fasta):
        self._outputDir = outputDir
        self._fasta = fasta
        self._images = []

    def addImage(self, imageBasename, title, hitInfo):
        self._images.append({
            'hitInfo': hitInfo,
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
            hitInfo = image['hitInfo']
            readIds = self._writeFASTA(i, image)
            if len(hitInfo['features']) > 0:
                features = self._writeFeatures(i, image)
            fp.write("""
      <a id="small_%d"></a>
      <h3>%d: %s</h3>
      <p>
        <a href="#big_%d"><img src="%s" class="thumbnail"/></a>
        Number of reads that hit overall: %d.
        Number of HSPs in this selection: %d.
        <br/><a href="#big_%d">Full size image</a>.
"""
                     % (i, i, title, i, image['imageBasename'],
                        hitInfo['hitCount'], len(hitInfo['items']), i))

            url = NCBISequenceLinkURL(title)
            if url:
                fp.write("""\
        <br/><a href="%s" target="_blank">NCBI info on this target</a>.
"""
                         % url)
            if len(hitInfo['features']):

                fp.write("""\
        <br/><a href="%s">Features</a>
"""
                         % features)

            if len(hitInfo['items']):
                fp.write("""\
        <br/>Reads: <span class="reads">%s</span>
"""
                         % ', '.join(readIds))

            fp.write("""
        <br clear="all"/>
      </p>
""")

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
            for item in image['hitInfo']['items']:
                sequenceId = item['sequenceId']
                reads.append(self._fasta[sequenceId])
                ids.append(self._fasta[sequenceId].description)
            SeqIO.write(reads, fp, 'fasta')
        return ids

    def _writeFeatures(self, i, image):
        """
        If there are too many features, write a txt file containing
        the features as a table.

        i: The number of the image in self._images.
        image: A member of self._images.
        """
        filename = '%s/features-%d.txt' % (self._outputDir, i)
        featureList = image['hitInfo']['features']
        with open(filename, 'w') as fp:
            for item in featureList:
                items = str(item)
                fp.write(items + '\n')
        return filename

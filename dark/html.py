from IPython.display import HTML

from dark.reads import Reads


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


def _sortHTML(titlesAlignments, by, limit=None):
    """
    Return an C{IPython.display.HTML} object with the alignments sorted by the
    given attribute.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param by: A C{str}, one of 'length', 'maxScore', 'medianScore',
        'readCount', or 'title'.
    @param limit: An C{int} limit on the number of results to show.
    @return: An HTML instance with sorted titles and information about
        hit read count, length, and e-values.
    """
    out = []
    for i, title in enumerate(titlesAlignments.sortTitles(by), start=1):
        if limit is not None and i > limit:
            break
        titleAlignments = titlesAlignments[title]
        link = NCBISequenceLink(title, title)
        out.append(
            '%3d: reads=%d, len=%d, max=%s median=%s<br/>'
            '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%s' %
            (i, titleAlignments.readCount(), titleAlignments.subjectLength,
             titleAlignments.bestHsp().score.score,
             titleAlignments.medianScore(), link))
    return HTML('<pre>' + '<br/>'.join(out) + '</pre>')


def summarizeTitlesByTitle(titlesAlignments, limit=None):
    """
    Sort match titles by title

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param limit: An C{int} limit on the number of results to show.
    @return: An C{IPython.display.HTML} instance with match titles sorted by
        title.
    """
    return _sortHTML(titlesAlignments, 'title', limit)


def summarizeTitlesByCount(titlesAlignments, limit=None):
    """
    Sort match titles by read count.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param limit: An C{int} limit on the number of results to show.
    @return: An C{IPython.display.HTML} instance with match titles sorted by
        read count.
    """
    return _sortHTML(titlesAlignments, 'readCount', limit)


def summarizeTitlesByLength(titlesAlignments, limit=None):
    """
    Sort match titles by sequence length.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param limit: An C{int} limit on the number of results to show.
    @return: An C{IPython.display.HTML} instance with match titles sorted by
        sequence length.
    """
    return _sortHTML(titlesAlignments, 'length', limit)


def summarizeTitlesByMaxScore(titlesAlignments, limit=None):
    """
    Sort hit titles by maximum score.

    @param titlesAlignments: A L{dark.blast.BlastMatchs} instance.
    @param limit: An C{int} limit on the number of results to show.
    @return: An C{IPython.display.HTML} instance with hit titles sorted by
        max score.
    """
    return _sortHTML(titlesAlignments, 'maxScore', limit)


def summarizeTitlesByMedianScore(titlesAlignments, limit=None):
    """
    Sort match titles by median score.

    @param titlesAlignments: A L{dark.titles.TitlesAlignments} instance.
    @param limit: An C{int} limit on the number of results to show.
    @return: An C{IPython.display.HTML} instance with match titles sorted by
        median score.
    """
    return _sortHTML(titlesAlignments, 'medianScore', limit)


class AlignmentPanelHTML(object):
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

    def addImage(self, imageBasename, title, graphInfo):
        self._images.append({
            'graphInfo': graphInfo,
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
            titleAlignments = self._titlesAlignments[title]
            graphInfo = image['graphInfo']
            reads = self._writeFASTA(i, image)
            fp.write("""
      <a id="small_%d"></a>
      <h3>%d: %s</h3>
      <p>
        <a href="#big_%d"><img src="%s" class="thumbnail"/></a>
        Sequence length: %d base pairs.<br/>
        Number of reads that hit overall: %d.<br/>
        Number of HSPs: %d.<br/>
        <a href="#big_%d">Full size image</a>.
"""
                     % (i, i, title, i, image['imageBasename'],
                        titleAlignments.subjectLength,
                        titleAlignments.readCount(),
                        titleAlignments.hspCount(), i))

            url = NCBISequenceLinkURL(title)
            if url:
                fp.write('<br/><a href="%s" target="_blank">NCBI info on this '
                         'target</a>.' % url)

            # Write out feature information.
            if graphInfo['features'] is None:
                fp.write('<br/>Feature lookup was False (or we were offline).')
            elif len(graphInfo['features']) == 0:
                fp.write('<br/>There were no features.')
            else:
                fp.write('<br/><a href="%s">Features</a>' %
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

            if len(reads):
                fp.write('<br/>Reads: <span class="reads">%s</span>'
                         % ', '.join(read.id for read in reads))

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

        @param i: The number of the image in self._images.
        @param image: A member of self._images.
        @return: A C{dark.reads.Reads} instance holding the reads for the
            image title.
        """
        reads = Reads()
        title = image['title']
        titleAlignments = self._titlesAlignments[title]
        for titleAlignment in titleAlignments:
            reads.add(titleAlignment.read)
        filename = '%s/%d.fasta' % (self._outputDir, i)
        reads.save(filename, 'fasta')
        return reads

    def _writeFeatures(self, i, image):
        """
        Write a txt file containing the features as a table.

        @param i: The number of the image in self._images.
        @param image: A member of self._images.
        @return: The C{str} features file name - just the base name, not
            including the path to the file.
        """
        basename = 'features-%d.txt' % i
        filename = '%s/%s' % (self._outputDir, basename)
        featureList = image['graphInfo']['features']
        with open(filename, 'w') as fp:
            for feature in featureList:
                fp.write('%s\n\n' % feature.feature)
        return basename

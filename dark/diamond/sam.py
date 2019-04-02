from tempfile import TemporaryFile
from functools import partial

from dark import __version__ as VERSION
from dark.btop import btop2cigar
from dark.diamond.conversion import DiamondTabularFormat
from dark.proteins import SqliteIndex
from dark.reads import DNARead


class SingleStreamSAMWriter(object):
    """
    A class that collects SAM output lines and writes them to a single stream.

    @param ram: Use RAM (i.e., not a temporary file) to hold the non-header SAM
        output. This will run faster but use more memory since all non-header
        SAM output will be stored in RAM and only written out when the full
        header can be determined (when 'save' is called).
    @param progName: The C{str} name of the program that produced the SAM. This
        will appear in the SAM header.
    """
    def __init__(self, ram, progName):
        self._ram = ram
        self._progName = progName
        self._referenceNames = set()

        if ram:
            self._nonHeaderLines = []
            self._emit = self._nonHeaderLines.append
        else:
            self._tf = TemporaryFile(mode='w+t', encoding='utf-8')
            self._emit = partial(print, file=self._tf)

    def addMatch(self, line, referenceName):
        """
        Add a line of SAM output.

        @param line: A C{str} line of SAM output.
        @param referenceName: The C{str} name of the reference that this SAM
            line refers to (this could also be extracted from the SAM line).
        """
        self._emit(line)
        self._referenceNames.add(referenceName)

    def save(self, fp, referenceLengths):
        """
        Write SAM output to a file.

        @param fp: An open file pointer to write to.
        @param referenceLengths: A C{dict} mapping C{str} reference sequence
            ids to their C{int} lengths.
        """
        # Print SAM headers.
        print('\n'.join(
            [
                '@PG\tID:DIAMOND\tPN:DIAMOND',
                '@PG\tID:%s\tPN:%s (version %s)\tPP:DIAMOND' %
                (self._progName, self._progName, VERSION,),
                '@CO\t%s is from the dark-matter package '
                '(https://github.com/acorg/dark-matter/)' % self._progName,
            ] +
            [
                '@SQ\tSN:%s\tLN:%d' % (name, referenceLengths[name])
                for name in sorted(self._referenceNames)
            ]), file=fp)

        # Print non-header lines.
        if self._ram:
            print('\n'.join(self._nonHeaderLines), file=fp)
        else:
            tf = self._tf
            tf.seek(0)
            for line in tf:
                print(line, end='', file=fp)
            tf.close()


class DiamondSAMWriter(object):
    """
    A class that can convert DIAMOND tabular output to SAM.

    @param genomesProteinsDatabaseFilename: A ...
    @param perReferenceOutput: If C{True}, a separate SAM file is written for
        matches against each reference sequence. The files will have names
        corresponding to the accession number of the reference, with a .sam
        suffix and will be put into the directory specified by
        C{outputDirectory}.
    @param outputDirectory: The C{str} name of the directory to create
        per-reference SAM files in. This is only used when
        C{perReferenceOutput} is C{True}, otherwise the output is written to
        stdout as a single SAM file.
    @param mappingQuality: The mapping quality to use for the SAM MAPQ field
        (#5). The default (255) indicates that mapping quality information is
        not available.
    @param ram: Do not use a temporary file to hold the non-header SAM output.
        This will run faster but use more memory since all non-header SAM
        output will be stored in RAM and only written out when the full
        header can be determined.
    @param keepDescriptions: If C{True}, do not discard text after the first
        space in query or subject sequence ids. Note that this violates the
        SAM specification, but since SAM files are TAB-separated there may be
        only a small chance this will cause problems downstream.
    """
    def __init__(self, genomesProteinsDatabaseFilename,
                 perReferenceOutput=False, outputDirectory='.',
                 mappingQuality=255, ram=False, keepDescriptions=False):
        self._genomesProteinsDatabaseFilename = genomesProteinsDatabaseFilename
        self._genomesProteins = SqliteIndex(genomesProteinsDatabaseFilename)
        self._perReferenceOutput = perReferenceOutput
        self._outputDirectory = outputDirectory
        self._mappingQuality = mappingQuality
        self._ram = ram
        self._keepDescriptions = keepDescriptions
        self._strToDiamondDict = DiamondTabularFormat().diamondFieldsToDict
        self._referenceLengths = {}

        if ram:
            self._nonHeaderLines = []
            self._emit = self._nonHeaderLines.append
        else:
            tf = TemporaryFile(mode='w+t', encoding='utf-8')
            self._emit = partial(print, file=tf)

    def addMatch(self, diamondStr):
        """
        Add information from a row of DIAMOND tabular output.

        @param diamondStr: A C{str} with TAB-separated fileds from DIAMOND
            output format 6.
        """
        match = self._strToDiamondDict(diamondStr)

        idOnly = not self._keepDescriptions
        qseqid = match['qseqid'].split()[0] if idOnly else match['qseqid']
        stitle = match['stitle'].split()[0] if idOnly else match['stitle']

        genome = self._genomesProteins.findProtein(stitle)
        if not genome:
            raise ValueError(
                'Could not find subject %r in genomes and '
                'proteins database file %r.' %
                (stitle, self._genomesProteinsDatabaseFilename))

        self._referenceLengths[stitle] = genome['length']

        # If the query frame is less than zero, the match was with a
        # reverse complemented translation of the query. Put the reverse
        # complement into the SAM output, which seems to be standard /
        # accepted practice based on my web searches.  See e.g.,
        # https://www.biostars.org/p/131891/ for what Bowtie2 does and for
        # some comments on this issue for SAM/BAM files in general.
        if match['qframe'] > 0:
            flag = 0
            qseq = match['qseq']
            qqual = match['qqual'] or '*'
        else:
            flag = 16
            qseq = DNARead('id', match['qseq']).reverseComplement().sequence
            qqual = match['qqual'][::-1] if match['qqual'] else '*'

        # Make a CIGAR string, including hard-clipped bases at the start and
        # end of the query (DIAMOND outputs a hard-clipped query sequence).
        startClipCount = match['qstart'] - 1
        endClipCount = match['qlen'] - match['qend']

        assert startClipCount >= 0
        assert endClipCount >= 0, (
            'Query sequence %s has length %d but the qend value is %d' %
            (qseq, len(match['qseq']), match['qend']))

        cigar = (
            ('%dH' % startClipCount if startClipCount else '') +
            btop2cigar(match['btop'], concise=False, aa=True) +
            ('%dH' % endClipCount if endClipCount else ''))

        self._emit('\t'.join(map(str, [
            # 1. QNAME
            qseqid,
            # 2. FLAG
            flag,
            # 3. RNAME
            stitle,
            # 4. POS. This needs to be a 1-based offset into the
            # nucleotide-equivalent of the DIAMOND subject sequence (which was
            # a protein since that is how DIAMOND operates). Because DIAMOND
            # gives back a 1-based protein location, we adjust to 0-based,
            # multiply by 3 to get to nucleotides, then adjust to 1-based.
            3 * (match['sstart'] - 1) + 1,
            # 5. MAPQ
            self._mappingQuality,
            # 6. CIGAR
            cigar,
            # 7. RNEXT
            '*',
            # 8. PNEXT
            0,
            # 9. TLEN
            0,
            # 10. SEQ
            qseq,
            # 11. QUAL
            qqual,
            # 12. Alignment score
            'AS:i:%d' % int(match['bitscore'])])))

    def save(self, fp):
        """
        Write SAM output to a file.
        """

from __future__ import print_function

from os import unlink
from tempfile import TemporaryFile
from resource import getrlimit, RLIMIT_NOFILE

from dark.btop import btop2cigar
from dark.diamond.conversion import DiamondTabularFormat
from dark.fpcache import FilePointerCache
from dark.genbank import GenomeRanges
from dark.reads import DNARead
from dark.utils import asHandle


class SAMWriter(object):
    """
    A class that collects SAM output lines and writes SAM to a file.

    @param ram: If C{True}, do not use a temporary file to hold the non-header
        SAM output. This will run faster but use more memory since all
        non-header SAM output will be stored in RAM and only written out when
        the full header can be determined.
    """
    def __init__(self, ram=False):
        self._ram = ram
        self._referenceLengths = {}

    def addMatchLine(self, line, genome, fp):
        """
        Add a line of SAM output.

        @param line: A C{str} line of valid pre-formatted SAM output.
        @param genome: A C{dict} with information about the (nucleotide) genome
            that the protein in the DIAMOND match comes from. The C{dict} is
            as returned by C{dark.proteins.SqliteIndex.findGenome}.
        @param fp: The file pointer to write the line to.
        """
        genomeName = genome['accession']
        genomeLength = genome['length']
        try:
            preexistingLength = self._referenceLengths[genomeName]
        except KeyError:
            self._referenceLengths[genomeName] = genomeLength
        else:
            if preexistingLength != genomeLength:
                raise ValueError(
                    'Reference %r passed with length %d but has already been '
                    'given with a different length (%d).' %
                    (genomeName, genomeLength, preexistingLength))

        print(line, file=fp)

    def save(self, fp, tmpfp=None, headerLines=None):
        """
        Write SAM output to a file.

        @param fp: An open file pointer to write the full SAM to.
        @param tmpfp: An open file pointer with the non-header SAM output.
            Only used if C{self._ram} is C{False}
        @param headerLines: An optional iterable of C{str}s SAM header
            lines (with no trailing newlines).
        """
        # Print SAM headers.
        print('\n'.join(
            (headerLines or []) +
            [
                '@SQ\tSN:%s\tLN:%d' % (name, self._referenceLengths[name])
                for name in sorted(self._referenceLengths)
            ]), file=fp)

        # Print SAM match lines.
        if self._ram:
            # Trade speed for memory and print the non-header lines one at a
            # time instead of making and printing a potentially massive string
            # with '\n'.join() (although arguably we could just make the single
            # string, seeing as self._ram is True.)
            for line in self._nonHeaderLines:
                print(line, file=fp)
        else:
            tmpfp.seek(0)
            for line in tmpfp:
                print(line, end='', file=fp)


class _DiamondSAMWriter(object):
    """
    A base class for converting DIAMOND tabular output to SAM.

    @param genomesProteins: An open sqlite3 database with information on
        genomes and proteins, as made by the dark-matter
        make-protein-database.py script.
    @param mappingQuality: The mapping quality to use for the SAM MAPQ field
        (#5). The default (255) indicates that mapping quality information is
        not available.
    @param ram: If C{True}, do not use a temporary file to hold the non-header
        SAM output. This will run faster but use more memory since all
        non-header SAM output will be stored in RAM and only written out when
        the full header can be determined.
    @param keepDescriptions: If C{True}, do not discard text after the first
        space in query or subject sequence ids. Note that this violates the
        SAM specification, but since SAM files are TAB-separated there may be
        only a small chance this will cause problems downstream.
    """
    def __init__(self, genomesProteins, mappingQuality=255, ram=False,
                 keepDescriptions=False):
        self._genomesProteins = genomesProteins
        self._mappingQuality = mappingQuality
        self._ram = ram
        self._keepDescriptions = keepDescriptions
        self._strToDiamondDict = DiamondTabularFormat().diamondFieldsToDict
        self._referenceLengths = {}

    def _preprocessMatch(self, diamondStr):
        """
        Add information from a row of DIAMOND tabular output.

        @param diamondStr: A C{str} with TAB-separated fields from DIAMOND
            output format 6.
        @return: A C{tuple} containing the match C{dict} (from
            C{self._strToDiamondDict}), the protein C{dict} as looked up
            by C{self._genomesProteins.findProtein}, and the genome C{dict}
            as looked up by C{self._genomesProteins.findGenome}.
        """
        match = self._strToDiamondDict(diamondStr)
        stitle = match['stitle']

        genomeAccession = self._genomesProteins.genomeAccession(stitle)
        genome = self._genomesProteins.findGenome(genomeAccession)
        if genome is None:
            raise ValueError(
                'Could not find accession %r in genomes database table.' %
                genomeAccession)

        proteinAccession = self._genomesProteins.proteinAccession(stitle)
        protein = self._genomesProteins.findProtein(proteinAccession)
        if protein is None:
            raise ValueError(
                'Could not find accession %r in proteins database table.' %
                proteinAccession)

        self._referenceLengths[genomeAccession] = genome['length']

        return match, protein, genome

    def _SAMLine(self, match, protein, genome):
        """
        Convert DIAMOND match information to a line of SAM file output.

        @param match: A C{dict} with information about the DIAMOND match, as
            returned by C{DiamondTabularFormat().diamondFieldsToDict} which
            has been called for us by C{self._preprocessMatch}.
        @param protein: A C{dict} with information about the protein the
            DIAMOND was for. The C{dict} is as returned by
            C{dark.proteins.SqliteIndex.findProtein}.
        @param genome: A C{dict} with information about the (nucleotide) genome
            that the protein in the DIAMOND match comes from. The C{dict} is
            as returned by C{dark.proteins.SqliteIndex.findGenome}.
        @return: A TAB-separated C{str} line of SAM.
        """
        qseqid = (match['qseqid'] if self._keepDescriptions else
                  match['qseqid'].split(None, 1)[0])

        ranges = GenomeRanges(protein['offsets'])

        # matchStartInGenome = ranges.startInGenome(match)
        # queryStartInGenome = matchStartInGenome - match['qstart'] - 1
        # orientations = ranges.orientations()

        print('ranges', ranges)
        print('protein', protein)
        print('match', match)
        # If the query frame is less than zero, the match was with a
        # translation of the reverse-complemented query. We'll put the
        # reverse complement into the SAM output. This seems to be standard
        # / accepted practice, based on my web searches.  See e.g.,
        # https://www.biostars.org/p/131891/ for what Bowtie2 does and for
        # some comments on this issue for SAM/BAM files in general.
        if match['qframe'] > 0:
            flag = 0
            qseq = match['full_qseq']
            qqual = match['full_qqual'] or '*'
        else:
            flag = 16
            qseq = DNARead('id',
                           match['full_qseq']).reverseComplement().sequence
            qqual = match['full_qqual'][::-1] if match['full_qqual'] else '*'

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

        return '\t'.join(map(str, (
            # 1. QNAME
            qseqid,
            # 2. FLAG
            flag,
            # 3. RNAME
            genome['accession'],
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
            'AS:i:%d' % int(match['bitscore']))))

    def addMatch(self, diamondStr):
        """
        Add information from a row of DIAMOND tabular output.

        @param diamondStr: A C{str} with TAB-separated fields from DIAMOND
            output format 6.
        """
        raise NotImplementedError('addMatch must be implemented by a subclass')

    def save(self, fp):
        """
        Write SAM output to a file.
        """
        raise NotImplementedError('save must be implemented by a subclass')


class SimpleDiamondSAMWriter(_DiamondSAMWriter):
    """
    Convert DIAMOND tabular output to SAM.

    @param genomesProteins: An open sqlite3 database with information on
        genomes and proteins, as made by the dark-matter
        make-protein-database.py script.
    @param mappingQuality: The mapping quality to use for the SAM MAPQ field
        (#5). The default (255) indicates that mapping quality information is
        not available.
    @param ram: If C{True}, do not use a temporary file to hold the non-header
        SAM output. This will run faster but use more memory since all
        non-header SAM output will be stored in RAM and only written out when
        the full header can be determined.
    @param keepDescriptions: If C{True}, do not discard text after the first
        space in query or subject sequence ids. Note that this violates the
        SAM specification, but since SAM files are TAB-separated there may be
        only a small chance this will cause problems downstream.
    """
    def __init__(self, genomesProteins, mappingQuality=255, ram=False,
                 keepDescriptions=False):
        _DiamondSAMWriter.__init__(
            self, genomesProteins, mappingQuality=mappingQuality, ram=ram,
            keepDescriptions=keepDescriptions)
        self._tf = TemporaryFile(mode='w+t', encoding='utf-8')
        self._writer = SAMWriter(ram=ram)

    def addMatch(self, diamondStr):
        """
        Add information from a row of DIAMOND tabular output.

        @param diamondStr: A C{str} with TAB-separated fields from DIAMOND
            output format 6.
        """
        match, protein, genome = self._preprocessMatch(diamondStr)
        self._writer.addMatchLine(
            self._SAMLine(match, protein, genome), genome, self._tf)

    def save(self, filename):
        """
        Write SAM output to a file.

        @param filename: A C{str} filename or an already-open file handle.
        """
        with asHandle(filename) as fp:
            self._writer.save(fp, tmpfp=self._tf)
        self._tf.close()
        self._tf = self._writer = None


class PerReferenceDiamondSAMWriter(_DiamondSAMWriter):
    """
    Convert DIAMOND tabular output to SAM and write it to per-reference
    output files.

    @param genomesProteins: An open sqlite3 database with information on
        genomes and proteins, as made by the dark-matter
        make-protein-database.py script.
    @param baseFilenameFunc: A function that produces a base filename to save
        SAM output to for a given reference name. A '.sam' suffix will be
        added to these files.
    @param mappingQuality: The mapping quality to use for the SAM MAPQ field
        (#5). The default (255) indicates that mapping quality information is
        not available.
    @param ram: If C{True}, do not use a temporary file to hold the non-header
        SAM output. This will run faster but use more memory since all
        non-header SAM output will be stored in RAM and only written out when
        the full header can be determined.
    @param keepDescriptions: If C{True}, do not discard text after the first
        space in query or subject sequence ids. Note that this violates the
        SAM specification, but since SAM files are TAB-separated there may be
        only a small chance this will cause problems downstream.
    @param fpcMaxsize: An C{int} maximum size for the file pointer cache used
        because SAM output may be being written to an arbitrary number of
        output files so we need to make sure the OS limit on open files is not
        exceeded and that we are not inefficiently constantly re-opening and
        writing single lines to files). If C{None} or C{0} is passed, a size
        of half the maximum number of files will be used.
    """
    SAM_SUFFIX = '.sam'
    TMP_SUFFIX = '.tmp'

    def __init__(self, genomesProteins, baseFilenameFunc=None,
                 mappingQuality=255, ram=False, keepDescriptions=False,
                 fpcMaxsize=0):
        _DiamondSAMWriter.__init__(
            self, genomesProteins, mappingQuality=mappingQuality, ram=ram,
            keepDescriptions=keepDescriptions)

        self._writers = {}
        self._baseFilenameFunc = baseFilenameFunc
        self._fpc = FilePointerCache(
            maxsize=(fpcMaxsize or getrlimit(RLIMIT_NOFILE)[0] >> 1),
            openArgs={'mode': 'wt'}, reopenArgs={'mode': 'at'})

    def addMatch(self, diamondStr):
        """
        Add information from a row of DIAMOND tabular output.

        @param diamondStr: A C{str} with TAB-separated fields from DIAMOND
            output format 6.
        """
        match, protein, genome = self._preprocessMatch(diamondStr)
        genomeAccession = genome['accession']
        basename = self._baseFilenameFunc(match['stitle'])

        try:
            writer, _ = self._writers[genomeAccession]
        except KeyError:
            writer = SAMWriter()
            self._writers[genomeAccession] = (writer, basename)

        writer.addMatchLine(self._SAMLine(match, protein, genome), genome,
                            self._fpc.open(basename + self.TMP_SUFFIX))

    def save(self):
        """
        Write SAM output to all files.
        """
        # Close the open file descriptor cache.
        self._fpc.close()

        # Write out SAM (including headers) for each reference.
        for writer, basename in self._writers.values():
            tmp = basename + self.TMP_SUFFIX
            with open(basename + self.SAM_SUFFIX, 'wt') as samfp:
                with open(tmp) as tmpfp:
                    writer.save(samfp, tmpfp)
            unlink(tmp)

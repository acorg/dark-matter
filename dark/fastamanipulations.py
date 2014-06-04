import sys
from Bio import SeqIO
from dark.conversion import JSONRecordsReader, XMLRecordsReader


def _getReadIds(readIdFile):
    """
    Given a file with readIds, outputs a list with all readIds, which can be
    read by fastaSubset().

    @param readIdFile: a C{str} of a filename with readIds.
    """
    readIds = []
    with open(readIdFile, 'r') as fp:
        for line in fp:
            readId = line.rstrip()
            readIds.append(readId)

    return readIds


def fastaSubset(inFile, readIds, present=True):
    """
    Given a set of reads and a fastafile, write reads to a new
    fastafile.

    @param inFile: a C{str} fastafile, from which reads should be subtracted.
    @param readIds: a C{list} of read identifiers
    @param present: a C{bool}. If True, the reads present in readIds will be
        written to outFile. If False, the reads NOT present in readIds will
        be written to outFile
    """
    # Test if readIds is a list, in which case it can be used directly, or if
    # a file has to be opened and transformed into a list.
    if type(readIds) == str:
        wanted = _getReadIds(readIds)
    else:
        wanted = readIds

    found = []
    for read in SeqIO.parse(inFile, 'fasta'):
        if present:
            if read.id in wanted:
                found.append(read)
                wanted.remove(read.id)
        else:
            if read.id not in wanted:
                found.append(read.id)
            if read.id in wanted:
                wanted.remove(read.id)

    if found:
        print >>sys.stderr, 'Found %d sequences.' % len(found)
        return found

    if wanted:
        print >>sys.stderr, 'WARNING: %d sequence%s not found: %s' % (
            len(wanted), '' if len(wanted) == 1 else 's were',
            ' '.join(wanted))
        sys.exit(1)


def fastaSubtract(fastaFiles):
    """
    Given a list of fastafiles, write either the reads in common to a new file.

    @param fastaFiles: a C{list} of fasta filenames of a fastafile
    """

    allReadIds = []
    individualReadIds = []
    for fastaFile in fastaFiles:
        for read in SeqIO.parse(fastaFile, 'fasta'):
            allReadIds.append(read.id)
            if read.id not in individualReadIds:
                individualReadIds.append(read.id)

    commonReadIds = []
    total = len(fastaFiles)
    for readId in individualReadIds:
        if allReadIds.count(readId) == total:
            commonReadIds.append(readId)

    fastaSubset(fastaFiles[0], commonReadIds, present=True)


def getReadsIdsFromBlast(blastFilenames, eCutoff=None, bitCutoff=None):
    """
    Given files with BLAST output, get the readIds which hit above an
    eCutoff or bitCutoff.

    @param blastFilenames: Either a single C{str} filename or a C{list} of
        C{str} file names containing BLAST output. Files can either be XML
        (-outfmt 5) BLAST output file or our smaller (possibly bzip2
        compressed) converted JSON equivalent produced by
        C{bin/convert-blast-xml-to-json.py} from a BLAST XML file.
    @param fastaFilename: the C{str} file name containing the sequences that
        were given to BLAST as queries.
    @param eCutoff: A C{float} e-value. Hits with e-value greater than or
            equal to this will be ignored.
    @param bitCutoff: A C{int} bitScore. Hits with a bitscore smaller than or
            equal to this will be ignored.
    """

    def records(blastFile):
        """
        Extract all BLAST records (up to C{self.limit}, if not C{None}).

        @return: A generator that yields BioPython C{Bio.Blast.Record.Blast}
            instances.
        """
        if blastFile.endswith('.xml'):
            reader = XMLRecordsReader(blastFile)
        elif (blastFile.endswith('.json') or
              blastFile.endswith('.json.bz2')):
            reader = JSONRecordsReader(blastFile)
        else:
            raise ValueError('Unknown BLAST file suffix for file %r.' %
                             blastFile)

        for record in reader.records():
            yield record

    readIds = []
    for blastFile in blastFilenames:
        for readNum, record in records(blastFile):
            for index, alignment in enumerate(record.alignments):
                if eCutoff:
                    if alignment.hsps[0].expect <= eCutoff:
                        readIds.append(record.query)
                elif bitCutoff:
                    if alignment.hsps[0].bits >= bitCutoff:
                        readIds.append(record.query)

    return readIds


def writeReadIds(readIds, outFile):
    """
    Write a list of readIds to a file specified by outFile.

    @param readIds: a C{list} containing the readIds to be written to the file.
    @param outFile: a C{str} name of the file to write to.
    """
    with open(outFile, 'w') as fp:
        fp.write('\n'.join(readIds))

    print >>sys.stderr, 'Wrote %d readIds to %s.' % (len(readIds), outFile)

import sys
from Bio import SeqIO
from dark.conversion import JSONRecordsReader, XMLRecordsReader


def _getReadIds(readIdFile):
    """
    Given a fasta file, outputs a list with all readIds, which can be
    read by fastaSubset().

    @param readIdFile: a C{str} of a fasta file.
    """
    readIds = []
    for record in SeqIO.parse(readIdFile, 'fasta'):
        readId = record.id
        readIds.append(readId)

    return readIds


def fastaSubset(inFile, readIds, invert=False):
    """
    Given a set of reads and a fastafile, return a list of sequences. If
    invert=False, return sequences specified in readIds, else return
    sequences not present in readIds.

    @param inFile: a C{str} fastafile, from which reads should be subtracted.
    @param readIds: a C{list} of read identifiers
    @param invert: a C{bool}. If True, the reads present in readIds will be
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
        if not invert:
            if read.id in wanted:
                found.append(read)
                wanted.remove(read.id)
        else:
            if read.id not in wanted:
                found.append(read)
            if read.id in wanted:
                wanted.remove(read.id)

    if found:
        print >>sys.stderr, 'Found %d sequences.' % len(found)
        return found

    if wanted:
        print >>sys.stderr, 'WARNING: %d sequence%s not found: %s' % (
            len(wanted), '' if len(wanted) == 1 else 's were',
            ' '.join(wanted))


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

    found = fastaSubset(fastaFiles[0], commonReadIds, invert=True)
    return found


def writeReadIds(readIds, outFile):
    """
    Write a list of readIds to a file specified by outFile.

    @param readIds: a C{list} containing the readIds to be written to the file.
    @param outFile: a C{str} name of the file to write to.
    """
    with open(outFile, 'w') as fp:
        fp.write('\n'.join(readIds))

    print >>sys.stderr, 'Wrote %d readIds to %s.' % (len(readIds), outFile)

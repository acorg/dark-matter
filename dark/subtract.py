from collections import defaultdict
from Bio import SeqIO


def fastaSubtract(fastaFiles):
    """
    Given a list of fasta files, either write the reads in common
    (invert=False) or the reads that are different (invert=True)
    to a new file.

    @param fastaFiles: a C{list} of fasta filenames of a fastafile
    @param invert: a C{bool}. If False, the reads present in readIds will be
        written to outFile. If True, the reads NOT present in readIds will
        be written to outFile.
    """
    allReads = defaultdict(int)
    firstFile = fastaFiles.pop(0)
    for read in SeqIO.parse(firstFile, 'fasta'):
        readId = read.id
        allReads[readId] += 1

    commonReads = []
    for fastaFile in fastaFiles:
        for read in SeqIO.parse(fastaFile, 'fasta'):
            readId = read.id
            if readId in allReads and allReads[readId] == 1:
                commonReads.append(read)
                allReads[readId] += 1

    return commonReads

from Bio import SeqIO
from hashlib import md5

from dark.reads import Read, Reads


def fastaToList(fastaFilename):
    return list(SeqIO.parse(fastaFilename, 'fasta'))


def dedupFasta(sequences):
    """
    sequences: an iterator producing Bio.Seq sequences.

    return: a generator of sequences with no duplicates.
    """
    seen = set()
    for s in sequences:
        thisSeq = str(s.seq)
        if thisSeq not in seen:
            seen.add(thisSeq)
            yield s


def dePrefixAndSuffixFasta(sequences):
    """
    sequences: an iterator producing Bio.Seq sequences.

    return: a generator of sequences with no duplicates and no fully contained
        subsequences.
    """
    sequences = sorted(sequences, key=lambda s: len(s.seq), reverse=True)
    seen = set()
    for s in sequences:
        thisSeq = str(s.seq)
        thisHash = md5(thisSeq).digest()
        if thisHash not in seen:
            # Add prefixes.
            newHash = md5()
            for nucl in thisSeq:
                newHash.update(nucl)
                seen.add(newHash.digest())
            # Add suffixes.
            for start in xrange(len(thisSeq) - 1):
                seen.add(md5(thisSeq[start + 1:]).digest())
            yield s
    # print 'seen contains %d hashes. total bytes = %d' % (len(seen),
    # len(seen) * len(md5('x').digest()))


def fastaSubtract(fastaFiles):
    """
    Given a list of open file descriptors, each with FASTA content,
    remove the reads found in the 2nd, 3rd, etc files from the first file
    in the list.

    @param fastaFiles: a C{list} of FASTA filenames.
    @raises IndexError: if passed an empty list.
    @return: An iterator producing C{Bio.SeqRecord} instances suitable for
        writing to a file using C{Bio.SeqIO.write}.

    """
    reads = {}
    firstFile = fastaFiles.pop(0)
    for seq in SeqIO.parse(firstFile, 'fasta'):
        reads[seq.id] = seq

    for fastaFile in fastaFiles:
        for seq in SeqIO.parse(fastaFile, 'fasta'):
            # Make sure that reads with the same id have the same sequence.
            if seq.id in reads:
                assert str(seq.seq) == str(reads[seq.id].seq)
            reads.pop(seq.id, None)

    return reads.itervalues()


class FastaReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to FASTA reads.

    @param file_: A C{str} file name or file handle, containing
        sequences in FASTA format,
    """
    def __init__(self, file_):
        self.file_ = file_
        Reads.__init__(self)

    def iter(self):
        """
        Iterate over the sequences in self.file_, yielding each as a Read
        instance.
        """
        for seq in SeqIO.parse(self.file_, 'fasta'):
            yield Read(seq.id, str(seq.seq))

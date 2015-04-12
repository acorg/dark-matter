from Bio import SeqIO
from hashlib import md5

from dark.reads import Reads, DNARead


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
    @param readClass: The class of read that should be yielded by iter.
    """
    def __init__(self, file_, readClass=DNARead):
        self.file_ = file_
        self.readClass = readClass
        Reads.__init__(self)

    def iter(self):
        """
        Iterate over the sequences in self.file_, yielding each as an
        instance of the desired read class.
        """
        for seq in SeqIO.parse(self.file_, 'fasta'):
            read = self.readClass(seq.description, str(seq.seq))
            read.checkAlphabet()
            yield read


def combineReads(filename, sequences, readClass=DNARead,
                 idPrefix='command-line-read-'):
    """
    Combine FASTA reads from a file and/or sequence strings.

    @param filename: A C{str} file name containing FASTA reads.
    @param sequences: A C{list} of C{str} sequences. If a sequence
        contains spaces, the last field (after splitting on spaces) will be
        used as the sequence and the first fields will be used as the sequence
        id.
    @param readClass: The class of the individual reads.
    @param idPrefix: The C{str} prefix that will be used for the id of the
        sequences in C{sequences} that do not have an id specified. A trailing
        sequence number will be appended to this prefix. Note that
        'command-line-read-', the default id prefix, could collide with ids in
        the FASTA file, if given. So output might be ambiguous. That's why we
        allow the caller to specify a custom prefix.
    @return: A C{FastaReads} instance.
    """
    # Read sequences from a FASTA file, if given.
    if filename:
        reads = FastaReads(filename, readClass=readClass)
    else:
        reads = Reads()

    # Add any individually specified subject sequences.
    if sequences:
        for count, sequence in enumerate(sequences, start=1):
            # Try splitting the sequence on its last space and using the
            # first part of the split as the read id. If there's no space,
            # assign a generic id.
            parts = sequence.rsplit(' ', 1)
            if len(parts) == 2:
                readId, sequence = parts
            else:
                readId = '%s%d' % (idPrefix, count)
            read = readClass(readId, sequence)
            read.checkAlphabet()
            reads.add(read)

    return reads

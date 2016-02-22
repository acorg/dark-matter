import six
from hashlib import md5

from Bio import SeqIO

from dark.reads import Reads, DNARead


def fastaToList(fastaFilename):
    return list(SeqIO.parse(fastaFilename, 'fasta'))


def dedupFasta(reads):
    """
    Remove sequence duplicates (based on sequence) from FASTA.

    @param reads: a C{dark.reads.Reads} instance.
    @return: a generator of C{dark.reads.Read} instances with no duplicates.
    """
    seen = set()
    add = seen.add
    for read in reads:
        hash_ = md5(read.sequence.encode('UTF-8')).digest()
        if hash_ not in seen:
            add(hash_)
            yield read


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
        thisHash = md5(thisSeq.encode('UTF-8')).digest()
        if thisHash not in seen:
            # Add prefixes.
            newHash = md5()
            for nucl in thisSeq:
                newHash.update(nucl.encode('UTF-8'))
                seen.add(newHash.digest())
            # Add suffixes.
            for start in range(len(thisSeq) - 1):
                seen.add(md5(thisSeq[start + 1:].encode('UTF-8')).digest())
            yield s


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

    return iter(reads.values())


class FastaReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to FASTA reads.

    @param _file: A C{str} file name or file handle, containing
        sequences in FASTA format,
    @param readClass: The class of read that should be yielded by iter.
    @param checkAlphabet: An C{int} or C{None}. If C{None}, alphabet checking
        will be done on all reads. If an C{int}, only that many reads will be
        checked. (Pass zero to have no checks done.)
    @param upperCase: If C{True}, reads will be converted to upper case.
    """
    def __init__(self, _file, readClass=DNARead, checkAlphabet=None,
                 upperCase=False):
        self._file = _file
        self._readClass = readClass
        self._checkAlphabet = checkAlphabet
        # TODO: It would be better if upperCase were an argument that could
        # be passed to Reads.__init__ and that could do the uppercasing in
        # its add method (as opposed to using it below in our iter method).
        # In that case, in the iter of this class we'd call self.add on
        # each of the sequences coming from self._file. Or, if we'd already
        # read the file we'd return Reads.iter(self) to re-iterate over the
        # sequences already added from the file.
        self._upperCase = upperCase
        if six.PY3:
            super().__init__()
        else:
            Reads.__init__(self)

    def iter(self):
        """
        Iterate over the sequences in self.file_, yielding each as an
        instance of the desired read class.
        """
        checkAlphabet = self._checkAlphabet
        for count, seq in enumerate(SeqIO.parse(self._file, 'fasta')):
            if self._upperCase:
                read = self._readClass(seq.description,
                                       str(seq.seq.upper()))
            else:
                read = self._readClass(seq.description, str(seq.seq))
            if checkAlphabet is None or count < checkAlphabet:
                read.checkAlphabet(count=None)
            yield read


def combineReads(filename, sequences, readClass=DNARead,
                 upperCase=False, idPrefix='command-line-read-'):
    """
    Combine FASTA reads from a file and/or sequence strings.

    @param filename: A C{str} file name containing FASTA reads.
    @param sequences: A C{list} of C{str} sequences. If a sequence
        contains spaces, the last field (after splitting on spaces) will be
        used as the sequence and the first fields will be used as the sequence
        id.
    @param readClass: The class of the individual reads.
    @param upperCase: If C{True}, reads will be converted to upper case.
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
        reads = FastaReads(filename, readClass=readClass, upperCase=upperCase)
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
            if upperCase:
                sequence = sequence.upper()
            read = readClass(readId, sequence)
            read.checkAlphabet()
            reads.add(read)

    return reads

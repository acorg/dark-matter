import six
from Bio import SeqIO

from dark.reads import Reads, SSAARead


class SSFastaReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to reads in FASTA format
    along with their PDB secondary structure.

    The pdb file is in FASTA format. For each pair of sequence / structure,
    it gives the amino acid sequence in the first record and the predicted
    secondary structure in the next. This is the file format used to hold
    results from the DSSP program. Secondary structures found by DSSP can be
    downloaded from http://www.rcsb.org/pdb/files/ss.txt on 11/11/2015

    IMPORTANT NOTE: the ss.txt file contains spaces in the structure
    records.  SeqIO.parse will silently collapse these to nothing, which
    will result in unequal length sequence and structure strings. So you
    will need to replace the spaces in that file with something else, like
    '-', to make sure the structure information has the correct length and
    alignment with the sequence.

    @param _file: A C{str} file name or file handle, containing
        sequences and structures in FASTA format,
    @param readClass: The class of read that should be yielded by iter. This
        must accept 3 C{str} arguments: an id, the sequence, the structure.
    @param checkAlphabet: An C{int} or C{None}. If C{None}, alphabet checking
        will be done on all reads. If an C{int}, only that many reads will be
        checked. (Pass zero to have no checks done.)
    @param upperCase: If C{True}, both read and structure sequences will be
        converted to upper case.
    """
    def __init__(self, _file, readClass=SSAARead, checkAlphabet=None,
                 upperCase=False):
        self._file = _file
        self._readClass = readClass
        self._checkAlphabet = checkAlphabet
        self._upperCase = upperCase
        if six.PY3:
            super().__init__()
        else:
            Reads.__init__(self)

    def iter(self):
        """
        Iterate over the sequences in self.file_, yielding each as an
        instance of the desired read class.

        @raise ValueError: If the input file has an odd number of records or
            if any sequence has a different length than its predicted
            secondary structure.
        """
        checkAlphabet = self._checkAlphabet
        upperCase = self._upperCase
        count = 0
        records = SeqIO.parse(self._file, 'fasta')
        while True:
            try:
                record = next(records)
            except StopIteration:
                break

            count += 1

            try:
                structureRecord = next(records)
            except StopIteration:
                raise ValueError('Structure file %r has an odd number of '
                                 'records.' % self._file)

            if len(structureRecord) != len(record):
                raise ValueError(
                    'Sequence %r length (%d) is not equal to structure %r '
                    'length (%d) in input file %r.' % (
                        record.description, len(record),
                        structureRecord.description, len(structureRecord),
                        self._file))

            if upperCase:
                read = self._readClass(record.description,
                                       str(record.seq.upper()),
                                       str(structureRecord.seq.upper()))
            else:
                read = self._readClass(record.description,
                                       str(record.seq),
                                       str(structureRecord.seq))

            if checkAlphabet is None or count < checkAlphabet:
                read.checkAlphabet(count=None)

            yield read

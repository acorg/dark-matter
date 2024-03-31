from Bio import SeqIO  # type: ignore

from dark.reads import Reads, SSAARead
from dark.utils import asHandle


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

    @param _files: Either a single C{str} file name or file handle, or a
        C{list} of C{str} file names and/or file handles. Each file or file
        handle must contain sequences in PDB FASTA format (see above).
    @param readClass: The class of read that should be yielded by iter. This
        must accept 3 C{str} arguments: an id, the sequence, the structure.
    @param upperCase: If C{True}, both read and structure sequences will be
        converted to upper case.
    """

    def __init__(self, _files, readClass=SSAARead, upperCase=False):
        self._files = _files if isinstance(_files, (list, tuple)) else [_files]
        self._readClass = readClass
        self._upperCase = upperCase
        super().__init__()

    def iter(self):
        """
        Iterate over the sequences in self.file_, yielding each as an
        instance of the desired read class.

        @raise ValueError: If the input file has an odd number of records or
            if any sequence has a different length than its predicted
            secondary structure.
        """
        upperCase = self._upperCase
        for _file in self._files:
            with asHandle(_file) as fp:
                records = SeqIO.parse(fp, "fasta")
                while True:
                    try:
                        record = next(records)
                    except StopIteration:
                        break
                    try:
                        structureRecord = next(records)
                    except StopIteration:
                        raise ValueError(
                            "Structure file %r has an odd number of records." % _file
                        )

                    if len(structureRecord) != len(record):
                        raise ValueError(
                            "Sequence %r length (%d) is not equal to "
                            "structure %r length (%d) in input file %r."
                            % (
                                record.description,
                                len(record),
                                structureRecord.description,
                                len(structureRecord),
                                _file,
                            )
                        )

                    if upperCase:
                        read = self._readClass(
                            record.description,
                            str(record.seq.upper()),
                            str(structureRecord.seq.upper()),
                        )
                    else:
                        read = self._readClass(
                            record.description,
                            str(record.seq),
                            str(structureRecord.seq),
                        )

                    yield read

from typing import Iterable

from prseq import FastaReader

from dark import File
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

    @param _files: Either a single C{str} file name or file handle, or a
        C{list} of C{str} file names and/or file handles. Each file or file
        handle must contain sequences in PDB FASTA format (see above).
    @param readClass: The class of read that should be yielded by iter. This
        must accept 3 C{str} arguments: an id, the sequence, the structure.
    @param upperCase: If C{True}, both read and structure sequences will be
        converted to upper case.
    """

    def __init__(
        self, _files: Iterable[File], readClass=SSAARead, upperCase: bool = False
    ):
        if isinstance(_files, tuple):
            self._files = _files
        elif isinstance(_files, list):
            self._files = tuple(_files)
        else:
            self._files = (_files,)

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
            assert isinstance(_file, File)
            with asHandle(_file) as fp:
                records = FastaReader(fp)
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

                    if len(structureRecord.sequence) != len(record.sequence):
                        raise ValueError(
                            "Sequence %r length (%d) is not equal to "
                            "structure %r length (%d) in input file %r."
                            % (
                                record.id,
                                len(record.sequence),
                                structureRecord.id,
                                len(structureRecord.sequence),
                                _file,
                            )
                        )

                    if upperCase:
                        read = self._readClass(
                            record.id,
                            str(record.sequence.upper()),
                            str(structureRecord.sequence.upper()),
                        )
                    else:
                        read = self._readClass(
                            record.id,
                            str(record.sequence),
                            str(structureRecord.sequence),
                        )

                    yield read

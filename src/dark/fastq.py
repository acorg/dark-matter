from prseq import FastqReader

from dark.reads import DNARead, Reads


class FastqReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to FASTQ reads.

    @param files: Either a single C{str} file name or file handle, or a
        C{list} of C{str} file names and/or file handles. Each file or file
        handle must contain sequences in FASTQ format.
    @param readClass: The class of read that should be yielded by iter.
    """

    def __init__(self, files, readClass=DNARead):
        self.files = files if isinstance(files, (list, tuple)) else [files]
        self.readClass = readClass
        super().__init__()

    def iter(self):
        """
        Iterate over the sequences in the files in self.files, yielding each
        as an instance of the desired read class.
        """
        for _file in self.files:
            for record in FastqReader(_file):
                yield self.readClass(record.id, record.sequence, record.quality)

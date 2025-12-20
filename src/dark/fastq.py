import sys
from collections.abc import Iterable

from prseq import FastqReader

from dark import File
from dark.reads import DNARead, Reads


class FastqReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to FASTQ reads.

    @param files: Either a single C{str} file name or file handle, or a
        C{list} of C{str} file names and/or file handles. Each file or file
        handle must contain sequences in FASTQ format.
    @param readClass: The class of read that should be yielded by iter.
    """

    def __init__(self, files: File | Iterable[File] | None = None, readClass=DNARead):
        if files is None:
            self.files = (sys.stdin.buffer,)
        elif isinstance(files, (list, tuple)):
            self.files = tuple(files)
        else:
             self.files = (files,)
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

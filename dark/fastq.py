from six import PY3

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from dark.reads import Reads, DNARead
from dark.utils import asHandle


class FastqReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to FASTQ reads.

    @param _files: Either a single C{str} file name or file handle, or a
        C{list} of C{str} file names and/or file handles. Each file or file
        handle must contain sequences in FASTQ format.
    @param readClass: The class of read that should be yielded by iter.
    """
    def __init__(self, _files, readClass=DNARead):
        self._files = _files if isinstance(_files, (list, tuple)) else [_files]
        self.readClass = readClass
        if PY3:
            super().__init__()
        else:
            Reads.__init__(self)

    def iter(self):
        """
        Iterate over the sequences in the files in self.files_, yielding each
        as an instance of the desired read class.
        """
        for _file in self._files:
            with asHandle(_file) as fp:
                # Use FastqGeneralIterator because it provides access to
                # the unconverted quality string (i.e., it doesn't try to
                # figure out the numeric quality values, which we don't
                # care about at this point).
                for sequenceId, sequence, quality in FastqGeneralIterator(fp):
                    yield self.readClass(sequenceId, sequence, quality)

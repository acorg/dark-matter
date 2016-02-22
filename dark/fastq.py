import six

from Bio.File import as_handle
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from dark.reads import Reads, DNARead


class FastqReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to FASTQ reads.

    @param file_: A C{str} file name or file handle, containing
        sequences in FASTQ format,
    @param readClass: The class of read that should be yielded by iter.
    """
    def __init__(self, file_, readClass=DNARead):
        self.file_ = file_
        self.readClass = readClass
        if six.PY3:
            super().__init__()
        else:
            Reads.__init__(self)

    def iter(self):
        """
        Iterate over the sequences in self.file_, yielding each as an
        instance of the desired read class.
        """
        # Use FastqGeneralIterator because it provides access to the
        # unconverted quality string (i.e., it doesn't try to figure out
        # the numeric quality values, which we don't care about at this
        # point).
        with as_handle(self.file_) as fp:
            for sequenceId, sequence, quality in FastqGeneralIterator(fp):
                yield self.readClass(sequenceId, sequence, quality)

# from collections import OrderedDict


class Read(object):
    """
    Hold information about a single read.

    @param id: A C{str} describing the read.
    @param sequence: A C{str} of sequence information (might be
        nucleotides or proteins).
    @param quality: An optional C{str} of phred quality scores. If not C{None},
        it must be the same length as C{sequence}.
    """
    def __init__(self, id, sequence, quality=None):
        if quality is not None and len(quality) != len(sequence):
            raise ValueError(
                'Invalid read: sequence length (%d) != quality length (%d)' %
                (len(sequence), len(quality)))
        self.id = id
        self.sequence = sequence.upper()
        self.quality = quality

    def __eq__(self, other):
        return (self.id == other.id and
                self.sequence == other.sequence and
                self.quality == other.quality)

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return len(self.sequence)


class Reads(object):
    """
    Maintain a collection of sequence reads.
    """

    def __init__(self):
        self.additionalReads = []
        self._length = 0

    def add(self, read):
        """
        Add a read to this collection of reads.

        @param read: A C{Read} instance.
        """
        self.additionalReads.append(read)
        self._length += 1

    def __iter__(self):
        # Reset self._length because __iter__ may be called more than once
        # and we don't want to increment _length each time we iterate.
        self._length = len(self.additionalReads)
        try:
            iter = self.iter()
        except AttributeError:
            pass
        else:
            for read in iter:
                self._length += 1
                yield read

        for read in self.additionalReads:
            yield read

    def __len__(self):
        # Note that len reflects the number of reads that are currently
        # known. If self.__iter__ has not been exhausted, we will not know
        # the true number of reads.
        return self._length

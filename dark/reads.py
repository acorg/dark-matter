from Bio.Data.IUPACData import (
    ambiguous_dna_complement, ambiguous_rna_complement)


def _makeComplementTable(complementData):
    """
    Make a sequence complement table.

    @param complementData: A C{dict} whose keys and values are strings of
        length one. A key, value pair indicates a substitution that should
        be performed during complementation.
    @return: A 256 character string that can be used as a translation table
        by the C{translate} method of a Python string.
    """
    table = range(256)
    for _from, to in complementData.iteritems():
        table[ord(_from[0])] = ord(to[0])
    return ''.join(map(chr, table))


_TRANSLATION_TABLE = {
    'dna': _makeComplementTable(ambiguous_dna_complement),
    'rna': _makeComplementTable(ambiguous_rna_complement),
}


class Read(object):
    """
    Hold information about a single read.

    @param id: A C{str} describing the read.
    @param sequence: A C{str} of sequence information (might be
        nucleotides or proteins).
    @param quality: An optional C{str} of phred quality scores. If not C{None},
        it must be the same length as C{sequence}.
    @param type: A C{str}, either 'dna', 'rna', or 'aa' to indicate the type
        of the sequence, either DNA, RNA, or amino acids.
    @raise ValueError: if the length of the quality string (if any) does not
        match the length of the sequence, or if the passed type is invalid.
    """
    def __init__(self, id, sequence, quality=None, type='dna'):
        if quality is not None and len(quality) != len(sequence):
            raise ValueError(
                'Invalid read: sequence length (%d) != quality length (%d)' %
                (len(sequence), len(quality)))
        if type not in ('aa', 'dna', 'rna'):
            raise ValueError('Unknown sequence type %r' % type)

        self.id = id
        self.sequence = sequence.upper()
        self.quality = quality
        self.type = type

    def __eq__(self, other):
        return (self.id == other.id and
                self.sequence == other.sequence and
                self.quality == other.quality)

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return len(self.sequence)

    def reverseComplement(self):
        """
        Reverse complement a nucleotide sequence.

        @raise ValueError: if an attempt is made to reverse complement an amino
            acid sequence.
        @return: The reverse complemented sequence as a C{Read} instance.
        """
        if self.type == 'aa':
            raise ValueError('Cannot reverse complement an amino acid '
                             'sequence')
        quality = None if self.quality is None else self.quality[::-1]
        sequence = self.sequence.translate(_TRANSLATION_TABLE[self.type])[::-1]
        return Read(self.id, sequence, quality)


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

        # Look for an 'iter' method in a possible subclass.
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
        # Note that len reflects the number of reads that are currently known.
        # If self.__iter__ has not been exhausted, we will not know the true
        # number of reads.
        return self._length

    def save(self, filename, format_='fasta'):
        """
        Write the reads to C{filename} in the requested format.

        @param filename: Either a C{str} file name to save into (the file will
            be overwritten) or an open file descriptor (e.g., sys.stdout).
        @param format_: A C{str} format to save as, either 'fasta' or 'fastq'.
        @raise ValueError: if C{format_} is 'fastq' and a read with no quality
            is present, or if an unknown format is requested.
        """
        format_ = format_.lower()
        if format_ == 'fasta':
            toString = lambda read: '>%s\n%s\n' % (read.id, read.sequence)
        elif format_ == 'fastq':
            for read in self:
                if read.quality is None:
                    raise ValueError('Read %r has no quality information' %
                                     read.id)
            toString = lambda read: '@%s\n%s\n+%s\n%s\n' % (
                read.id, read.sequence, read.id, read.quality)
        else:
            raise ValueError("Save format must be either 'fasta' or 'fastq'.")

        if isinstance(filename, basestring):
            with open(filename, 'w') as fp:
                for read in self:
                    fp.write(toString(read))
        else:
            # We have a file-like object.
            for read in self:
                filename.write(toString(read))

    def filter(self, minLength=None, maxLength=None):
        """
        Filter a set of reads to produce a matching subset.

        Note: there are many additional filtering options that could be added,
        e.g., filtering on read id (whitelist, blacklist, regex, etc), GC %,
        and quality.

        @param minLength: The minimum acceptable length.
        @param maxLength: The maximum acceptable length.
        """
        result = Reads()
        for read in self:
            readLen = len(read)
            if ((minLength is None or readLen >= minLength) and
                    (maxLength is None or readLen <= maxLength)):
                result.add(read)
        return result

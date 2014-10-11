from Bio.Seq import translate
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
                self.quality == other.quality and
                self.type == other.type)

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return len(self.sequence)

    def toString(self, format_):
        """
        Convert the read to a string format.

        @param format_: Either 'fasta' or 'fastq'.
        @raise ValueError: if C{format_} is 'fastq' and the read has no quality
            information, or if an unknown format is requested.
        @return: A C{str} representing the read in the requested format.
        """
        if format_ == 'fasta':
            return '>%s\n%s\n' % (self.id, self.sequence)
        elif format_ == 'fastq':
            if self.quality is None:
                raise ValueError('Read %r has no quality information' %
                                 self.id)
            else:
                return '@%s\n%s\n+%s\n%s\n' % (
                    self.id, self.sequence, self.id, self.quality)
        else:
            raise ValueError("Format must be either 'fasta' or 'fastq'.")

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

    def translations(self):
        """
        Yield all six translations of a nucleotide sequence.

        @raise ValueError: If the read is already amino acids.
        @return: A generator that produces six L{TranslatedRead} instances.
        """
        if self.type == 'aa':
            raise ValueError('Cannot translate an amino acid sequence')

        rc = self.reverseComplement().sequence
        for reverseComplemented in False, True:
            for frame in 0, 1, 2:
                seq = rc if reverseComplemented else self.sequence
                # Get the suffix of the sequence for translation. I.e.,
                # skip 0, 1, or 2 initial bases, depending on the frame.
                # Note that this makes a copy of the sequence, which we can
                # then safely append 'N' bases to to adjust its length to
                # be zero mod 3.
                suffix = seq[frame:]
                lengthMod3 = len(suffix) % 3
                if lengthMod3:
                    suffix += ('NN' if lengthMod3 == 1 else 'N')
                yield TranslatedRead(self, translate(suffix), frame,
                                     reverseComplemented)

    def aaToProperties(self):
        """
        Translate an amino acid sequence to properties.

        @raise ValueError: If sequence type is 'dna', 'rna' or 'properties'.
        @return: A generator that produces six L{PropertiesRead} instances.
        """
        # H = hydrophobic
        # P = polar
        # B = basic
        # A = acidic

        SIMPLE = {'A': 'H', 'V': 'H', 'M': 'H', 'L': 'H', 'I': 'H', 'P': 'H',
                  'W': 'H', 'F': 'H', 'K': 'B', 'R': 'B', 'H': 'B', 'Y': 'P',
                  'T': 'P', 'Q': 'P', 'G': 'P', 'S': 'P', 'C': 'P', 'D': 'P',
                  'E': 'A', 'D': 'A'}

        if self.type in ('dna', 'rna'):
            raise ValueError('Cannot convert nucleotide sequence to '
                             'properties.')
        if self.type == 'properties':
            raise ValueError('Sequence is a properties sequence already.')

        aaSeq = self.sequence
        properties = ''.join([SIMPLE[base] for base in aaSeq])
        yield PropertiesRead(self, properties)


class PropertiesRead(Read):
    """
    Hold information about one AA->properties translation of a Read.

    @param originalRead: The original AA L{Read} instance from which
        this translation was obtained.
    @param sequence: The C{str} properties translated sequence.
    """
    def __init__(self, originalRead, sequence):
        newId = '%s-properties' % originalRead.id

        Read.__init__(self, newId, sequence, type='properties')
        self.originalRead = originalRead

    def __eq__(self, other):
        return (Read.__eq__(self, other) and
                self.originalRead == other.originalRead)

    def __ne__(self, other):
        return not self == other


class TranslatedRead(Read):
    """
    Hold information about one DNA->AA translation of a Read.

    @param originalRead: The original DNA or RNA L{Read} instance from which
        this translation was obtained.
    @param sequence: The C{str} AA translated sequence.
    @param frame: The C{int} frame, either 0, 1, or 2.
    @param reverseComplemented: A C{bool}, C{True} if the original sequence
        must be reverse complemented to obtain this AA sequence.
    """
    def __init__(self, originalRead, sequence, frame,
                 reverseComplemented=False):
        if frame not in (0, 1, 2):
            raise ValueError('Frame must be 0, 1, or 2')
        newId = '%s-frame%d%s' % (originalRead.id, frame,
                                  'rc' if reverseComplemented else '')
        Read.__init__(self, newId, sequence, type='aa')
        self.originalRead = originalRead
        self.frame = frame
        self.reverseComplemented = reverseComplemented

    def __eq__(self, other):
        return (Read.__eq__(self, other) and
                self.frame == other.frame and
                self.reverseComplemented == other.reverseComplemented and
                self.originalRead == other.originalRead)

    def __ne__(self, other):
        return not self == other

    def maximumORFLength(self):
        """
        Return the length of the longest (possibly partial) ORF in a translated
        read. The ORF may originate or terminate outside the sequence, which is
        why the length is just a lower bound.
        """
        return max(len(orf) for orf in self.sequence.split('*'))


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

        if isinstance(filename, basestring):
            with open(filename, 'w') as fp:
                for read in self:
                    fp.write(read.toString(format_))
        else:
            # We have a file-like object.
            for read in self:
                filename.write(read.toString(format_))

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

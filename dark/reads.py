from Bio.Seq import translate
from Bio.Data.IUPACData import (
    ambiguous_dna_complement, ambiguous_rna_complement)

from dark.aa import PROPERTIES, NONE


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
    return ''.join(map(chr, table)), complementData


class Read(object):
    """
    Hold information about a single read.

    @param id: A C{str} describing the read.
    @param sequence: A C{str} of sequence information (might be
        nucleotides or proteins).
    @param quality: An optional C{str} of phred quality scores. If not C{None},
        it must be the same length as C{sequence}.
    @raise ValueError: if the length of the quality string (if any) does not
        match the length of the sequence.
    """
    def __init__(self, id, sequence, quality=None):
        if quality is not None and len(quality) != len(sequence):
            raise ValueError(
                'Invalid read: sequence length (%d) != quality length (%d)' %
                (len(sequence), len(quality)))

        self.id = id
        try:
            self.sequence = sequence.upper()
        except AttributeError:
            self.sequence = sequence
        self.quality = quality

    def __eq__(self, other):
        return (self.id == other.id and
                self.sequence == other.sequence and
                self.quality == other.quality)

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


class _NucleotideRead(Read):
    """
    Holds methods to work with nucleotide (DNA and RNA) sequences.
    """
    def translations(self):
        """
        Yield all six translations of a nucleotide sequence.

        @return: A generator that produces six L{TranslatedRead} instances.
        """
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

    def reverseComplement(self):
        """
        Reverse complement a nucleotide sequence.

        @return: The reverse complemented sequence as an instance of the
            current class.
        """
        quality = None if self.quality is None else self.quality[::-1]
        sequence = self.sequence.translate(self.COMPLEMENT_TABLE[0])[::-1]
        return self.__class__(self.id, sequence, quality)


class DNARead(_NucleotideRead):
    """
    Hold information and methods to work with DNA reads.
    """
    COMPLEMENT_TABLE = _makeComplementTable(ambiguous_dna_complement)


class RNARead(_NucleotideRead):
    """
    Hold information and methods to work with RNA reads.
    """
    COMPLEMENT_TABLE = _makeComplementTable(ambiguous_rna_complement)


class AARead(Read):
    """
    Hold information and methods to work with AA reads.
    """
    def properties(self):
        """
        Translate an amino acid sequence to properties.

        @return: A generator that produces an L{PropertiesRead} instance.
        """
        return [PROPERTIES[aa] if aa in PROPERTIES else NONE
                for aa in self.sequence]


class TranslatedRead(AARead):
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
        Read.__init__(self, newId, sequence)
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

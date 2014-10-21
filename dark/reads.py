from os import unlink
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
    return ''.join(map(chr, table))


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
        sequence = self.sequence.translate(self.COMPLEMENT_TABLE)[::-1]
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

        @return: A generator yielding properties for the residues in the
            current sequence.
        """
        return (PROPERTIES.get(aa, NONE) for aa in self.sequence)

    def ORFs(self):
        """
        Find all ORFs in our sequence.

        @return: A generator that yields AAReadORF instances that correspond
            to the ORFs found in the AA sequence.
        """
        ORFStart = None
        openLeft = True

        for index, residue in enumerate(self.sequence):
            if residue == 'M':
                # Start codon.
                openLeft = False
            elif residue == '*':
                # Stop codon. Yield an ORF, if it has non-zero length.
                if ORFStart is not None and index - ORFStart > 0:
                    # The ORF has non-zero length.
                    yield AAReadORF(self, ORFStart, index, openLeft, False)
                    ORFStart = None
            else:
                if ORFStart is None:
                    ORFStart = index

        # End of sequence. Yield the final ORF if there is one and it has
        # non-zero length.
        length = len(self.sequence)
        if ORFStart is not None and length - ORFStart > 0:
            yield AAReadORF(self, ORFStart, length, openLeft, True)


class AAReadORF(AARead):
    """
    Hold information about an ORF from an AA read.

    @param originalRead: The original L{AARead} instance in which this ORF
        occurs.
    @param start: The C{int} offset where the ORF starts in the original read.
    @param stop: The Python-style C{int} offset of the end of the ORF in the
        original read. The final index is not included in the ORF.
    @param openLeft: A C{bool}. If C{True}, the ORF potentially begins before
        the sequence given in C{sequence}. I.e., the ORF-detection code started
        to examine a read assuming it was already in an ORF. If C{False}, a
        start codon was found preceeding this ORF.
    @param openRight: A C{bool}. If C{True}, the ORF potentially ends after
        the sequence given in C{sequence}. I.e., the ORF-detection code
        was in an ORF when it encountered the end of a read (so no stop codon
        was found). If C{False}, a stop codon was found in the read after this
        ORF.
    """
    def __init__(self, originalRead, start, stop, openLeft, openRight):
        if start < 0:
            raise ValueError('start offset (%d) less than zero' % start)
        if stop > len(originalRead):
            raise ValueError('stop offset (%d) > original read length (%d)' %
                             (stop, len(originalRead)))
        if start > stop:
            raise ValueError('start offset (%d) greater than stop offset (%d)'
                             % (start, stop))
        newId = '%s-%s%d:%d%s' % (originalRead.id,
                                  '(' if openLeft else '[',
                                  start, stop,
                                  ')' if openRight else ']')
        AARead.__init__(self, newId, originalRead.sequence[start:stop])
        self.start = start
        self.stop = stop
        self.openLeft = openLeft
        self.openRight = openRight


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
        AARead.__init__(self, newId, sequence)
        self.frame = frame
        self.reverseComplemented = reverseComplemented

    def __eq__(self, other):
        return (AARead.__eq__(self, other) and
                self.frame == other.frame and
                self.reverseComplemented == other.reverseComplemented)

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
        @return: C{self} in case our caller wants to chain the result into
            another method call.
        """
        format_ = format_.lower()

        if isinstance(filename, basestring):
            try:
                with open(filename, 'w') as fp:
                    for read in self:
                        fp.write(read.toString(format_))
            except ValueError:
                unlink(filename)
                raise
        else:
            # We have a file-like object.
            for read in self:
                filename.write(read.toString(format_))
        return self

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

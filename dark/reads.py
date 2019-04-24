import sys
import six
from os import unlink
from functools import total_ordering
from collections import Counter
from hashlib import md5
from random import uniform

from Bio.Seq import translate
from Bio.Data.IUPACData import (
    ambiguous_dna_complement, ambiguous_rna_complement)

from dark.aa import AA_LETTERS, NAMES as AA_NAMES
from dark.filter import TitleFilter
from dark.aa import PROPERTIES, PROPERTY_DETAILS, NONE


def _makeComplementTable(complementData):
    """
    Make a sequence complement table.

    @param complementData: A C{dict} whose keys and values are strings of
        length one. A key, value pair indicates a substitution that should
        be performed during complementation.
    @return: A 256 character string that can be used as a translation table
        by the C{translate} method of a Python string.
    """
    table = list(range(256))
    for _from, to in complementData.items():
        table[ord(_from[0].lower())] = ord(to[0].lower())
        table[ord(_from[0].upper())] = ord(to[0].upper())
    return ''.join(map(chr, table))


@total_ordering
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
    ALPHABET = None

    def __init__(self, id, sequence, quality=None):
        if quality is not None and len(quality) != len(sequence):
            raise ValueError(
                'Invalid read: sequence length (%d) != quality length (%d)' %
                (len(sequence), len(quality)))

        self.id = id
        self.sequence = sequence
        self.quality = quality

    def __eq__(self, other):
        return (self.id == other.id and
                self.sequence == other.sequence and
                self.quality == other.quality)

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        return ((self.id, self.sequence, self.quality) <
                (other.id, other.sequence, other.quality))

    def __len__(self):
        return len(self.sequence)

    def __hash__(self):
        """
        Calculate a hash key for a read.

        @return: The C{int} hash key for the read.
        """
        if self.quality is None:
            return hash(md5(self.id.encode('UTF-8') + b'\0' +
                            self.sequence.encode('UTF-8')).digest())
        else:
            return hash(md5(self.id.encode('UTF-8') + b'\0' +
                            self.sequence.encode('UTF-8') + b'\0' +
                            self.quality.encode('UTF-8')).digest())

    def __getitem__(self, item):
        sequence = self.sequence[item]
        quality = None if self.quality is None else self.quality[item]
        return self.__class__(self.id, sequence, quality)

    def toString(self, format_):
        """
        Convert the read to a string format.

        @param format_: Either 'fasta', 'fastq' or 'fasta-ss'.
        @raise ValueError: if C{format_} is 'fastq' and the read has no quality
            information, if C{format_} is 'fasta-ss' and the read has no
            structure information, or if an unknown format is requested.
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
            raise ValueError("Format must be either 'fasta', 'fastq' or "
                             "'fasta-ss'.")

    def toDict(self):
        """
        Get information about this read in a dictionary.

        @return: A C{dict} with keys/values for the attributes of self.
        """
        return {
            'id': self.id,
            'sequence': self.sequence,
            'quality': self.quality,
        }

    def reverse(self):
        """
        Reverse a read (note that this is NOT a reverse complement).

        @return: The reversed sequence as an instance of the current class.
        """
        return self.__class__(
            self.id,
            self.sequence[::-1],
            None if self.quality is None else self.quality[::-1])

    @classmethod
    def fromDict(cls, d):
        """
        Create a new instance from attribute values provided in a dictionary.

        @param d: A C{dict} with keys/values for the attributes of a new
            instance of this class. Keys 'id' and 'sequence' with C{str} values
            must be provided. A 'quality' C{str} key is optional.
        @return: A new instance of this class, with values taken from C{d}.
        """
        return cls(d['id'], d['sequence'], d.get('quality'))

    def lowComplexityFraction(self):
        """
        What fraction of a read's bases are in low-complexity regions?
        By convention, a region of low complexity is indicated by lowercase
        base letters.

        @return: The C{float} representing the fraction of bases in the
            read that are in regions of low complexity.
        """
        length = len(self)
        if length:
            lowerCount = len(list(filter(str.islower, self.sequence)))
            return float(lowerCount) / length
        else:
            return 0.0

    def walkHSP(self, hsp, includeWhiskers=True):
        """
        Provide information about exactly how a read matches a subject, as
        specified by C{hsp}.

        @param includeWhiskers: If C{True} yield information from the
            (possibly empty) non-matching ends of the read.
        @return: A generator that yields (offset, residue, inMatch) tuples.
            The offset is the offset into the matched subject. The residue is
            the base in the read (which might be '-' to indicate a gap in the
            read was aligned with the subject at this offset). inMatch will be
            C{True} for residues that are part of the HSP match, and C{False}
            for the (possibly non-existent) parts of the read that fall outside
            the HSP (aka, the "whiskers" in an alignment graph).
        """

        # It will be easier to understand the following implementation if
        # you refer to the ASCII art illustration of an HSP in the
        # dark.hsp._Base class in hsp.py

        # Left whisker.
        if includeWhiskers:
            readOffset = 0
            subjectOffset = hsp.readStartInSubject
            while subjectOffset < hsp.subjectStart:
                yield (subjectOffset, self.sequence[readOffset], False)
                readOffset += 1
                subjectOffset += 1

        # Match.
        for subjectOffset, residue in enumerate(hsp.readMatchedSequence,
                                                start=hsp.subjectStart):
            yield (subjectOffset, residue, True)

        # Right whisker.
        if includeWhiskers:
            readOffset = hsp.readEnd
            subjectOffset = hsp.subjectEnd
            while subjectOffset < hsp.readEndInSubject:
                yield (subjectOffset, self.sequence[readOffset], False)
                readOffset += 1
                subjectOffset += 1

    def checkAlphabet(self, count=10):
        """
        A function which checks whether the sequence in a L{dark.Read} object
        corresponds to its readClass. For AA reads, more testing is done in
        dark.Read.AARead.checkAlphabet.

        @param count: An C{int}, indicating how many bases or amino acids at
            the start of the sequence should be considered. If C{None}, all
            bases are checked.
        @return: C{True} if the alphabet characters in the first C{count}
            positions of sequence is a subset of the allowed alphabet for this
            read class, or if the read class has a C{None} alphabet.
        @raise ValueError: If the sequence alphabet is not a subset of the read
            class alphabet.
        """
        if count is None:
            readLetters = set(self.sequence.upper())
        else:
            readLetters = set(self.sequence.upper()[:count])
        # Check if readLetters is a subset of self.ALPHABET.
        if self.ALPHABET is None or readLetters.issubset(self.ALPHABET):
            return readLetters
        raise ValueError("Read alphabet (%r) is not a subset of expected "
                         "alphabet (%r) for read class %s." % (
                             ''.join(sorted(readLetters)),
                             ''.join(sorted(self.ALPHABET)),
                             str(self.__class__.__name__)))

    def newFromSites(self, sites, exclude=False):
        """
        Create a new read from self, with only certain sites.

        @param sites: A set of C{int} 0-based sites (i.e., indices) in
            sequences that should be kept. If C{None} (the default), all sites
            are kept.
        @param exclude: If C{True} the C{sites} will be excluded, not
            included.
        """
        if exclude:
            sites = set(range(len(self))) - sites

        newSequence = []
        if self.quality:
            newQuality = []
            for index, (base, quality) in enumerate(zip(self.sequence,
                                                        self.quality)):
                if index in sites:
                    newSequence.append(base)
                    newQuality.append(quality)
            read = self.__class__(self.id, ''.join(newSequence),
                                  ''.join(newQuality))
        else:
            for index, base in enumerate(self.sequence):
                if index in sites:
                    newSequence.append(base)
            read = self.__class__(self.id, ''.join(newSequence))

        return read


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
    ALPHABET = set('ATCG')

    COMPLEMENT_TABLE = _makeComplementTable(ambiguous_dna_complement)


class RNARead(_NucleotideRead):
    """
    Hold information and methods to work with RNA reads.
    """
    ALPHABET = set('ATCGU')

    COMPLEMENT_TABLE = _makeComplementTable(ambiguous_rna_complement)


class DNAKozakRead(DNARead):
    """
    Hold information about a Kozak sequence.

    @param originalRead: The C{dark.reads.Read} instance in which the Kozak
        sequence was found.
    @param start: The C{int} start location of the Kozak sequence.
    @param stop: The C{int} stop location of the Kozak sequence (this is a
        Python string index, so the final Kozak sequence character is the one
        before this offset in the sequence.
    @param kozakQuality: A C{float}, giving the percentage of the 5 variable
        locations in the Kozak sequence that match the most frequent Kozak
        nucleotides.
    """
    def __init__(self, originalRead, start, stop, kozakQuality):
        if start < 0:
            raise ValueError('start offset (%d) less than zero' % start)
        if stop > len(originalRead):
            raise ValueError('stop offset (%d) > original read length (%d)' %
                             (stop, len(originalRead)))
        if start > stop:
            raise ValueError('start offset (%d) greater than stop offset (%d)'
                             % (start, stop))

        newId = '%s-(%d:%d)' % (originalRead.id, start, stop)

        if originalRead.quality:
            DNARead.__init__(self, newId, originalRead.sequence[start:stop],
                             originalRead.quality[start:stop])
        else:
            DNARead.__init__(self, newId, originalRead.sequence[start:stop])
        self.originalRead = originalRead
        self.start = start
        self.stop = stop
        self.kozakQuality = kozakQuality

    def __eq__(self, other):
        return (self.id == other.id and
                self.sequence == other.sequence and
                self.originalRead == other.originalRead and
                self.quality == other.quality and
                self.start == other.start and
                self.stop == other.stop and
                self.kozakQuality == other.kozakQuality)


# Keep a single GOR4 instance that can be used by all AA reads. This saves us
# from re-scanning the GOR IV secondary structure database every time we make
# an AARead instance. This will be initialized when it's first needed.
_GOR4 = None


class AARead(Read):
    """
    Hold information and methods to work with AA reads.
    """
    ALPHABET = set(AA_LETTERS)

    def checkAlphabet(self, count=10):
        """
        A function which checks if an AA read really contains amino acids. This
        additional testing is needed, because the letters in the DNA alphabet
        are also in the AA alphabet.

        @param count: An C{int}, indicating how many bases or amino acids at
            the start of the sequence should be considered. If C{None}, all
            bases are checked.
        @return: C{True} if the alphabet characters in the first C{count}
            positions of sequence is a subset of the allowed alphabet for this
            read class, or if the read class has a C{None} alphabet.
        @raise ValueError: If a DNA sequence has been passed to AARead().
        """
        if six.PY3:
            readLetters = super().checkAlphabet(count)
        else:
            readLetters = Read.checkAlphabet(self, count)
        if len(self) > 10 and readLetters.issubset(set('ACGT')):
            raise ValueError('It looks like a DNA sequence has been passed to '
                             'AARead().')
        return readLetters

    def properties(self):
        """
        Translate an amino acid sequence to properties of the form:
        'F': HYDROPHOBIC | AROMATIC.

        @return: A generator yielding properties for the residues in the
            current sequence.
        """
        return (PROPERTIES.get(aa, NONE) for aa in self.sequence)

    def propertyDetails(self):
        """
        Translate an amino acid sequence to properties. Each property of the
        amino acid gets a value scaled from -1 to 1.

        @return: A list of property dictionaries.
        """
        return (PROPERTY_DETAILS.get(aa, NONE) for aa in self.sequence)

    def ORFs(self, openORFs=False):
        """
        Find all ORFs in our sequence.

        @param openORFs: If C{True} allow ORFs that do not have a start codon
            and/or do not have a stop codon.
        @return: A generator that yields AAReadORF instances that correspond
            to the ORFs found in the AA sequence.
        """

        # Return open ORFs to the left and right and closed ORFs within the
        # sequence.
        if openORFs:
            ORFStart = 0
            inOpenORF = True  # open on the left
            inORF = False

            for index, residue in enumerate(self.sequence):
                if residue == '*':
                    if inOpenORF:
                        if index:
                            yield AAReadORF(self, ORFStart, index, True, False)
                        inOpenORF = False
                    elif inORF:
                        if ORFStart != index:
                            yield AAReadORF(self, ORFStart, index,
                                            False, False)
                        inORF = False
                elif residue == 'M':
                    if not inOpenORF and not inORF:
                        ORFStart = index + 1
                        inORF = True

            # End of sequence. Yield the final ORF, open to the right, if
            # there is one and it has non-zero length.
            length = len(self.sequence)
            if inOpenORF and length > 0:
                yield AAReadORF(self, ORFStart, length, True, True)
            elif inORF and ORFStart < length:
                yield AAReadORF(self, ORFStart, length, False, True)

        # Return only closed ORFs.
        else:
            inORF = False

            for index, residue in enumerate(self.sequence):
                if residue == 'M':
                    if not inORF:
                        inORF = True
                        ORFStart = index + 1
                elif residue == '*':
                    if inORF:
                        if not ORFStart == index:
                            yield AAReadORF(self, ORFStart,
                                            index, False, False)
                        inORF = False


class AAReadWithX(AARead):
    """
    Hold information and methods to work with AA reads with additional
    characters.
    """
    ALPHABET = set(AA_LETTERS + ['X'])


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
        if originalRead.quality:
            AARead.__init__(self, newId, originalRead.sequence[start:stop],
                            originalRead.quality[start:stop])
        else:
            AARead.__init__(self, newId, originalRead.sequence[start:stop])
        self.start = start
        self.stop = stop
        self.openLeft = openLeft
        self.openRight = openRight

    def toDict(self):
        """
        Get information about this read in a dictionary.

        @return: A C{dict} with keys/values for the attributes of self.
        """
        if six.PY3:
            result = super().toDict()
        else:
            result = AARead.toDict(self)

        result.update({
            'start': self.start,
            'stop': self.stop,
            'openLeft': self.openLeft,
            'openRight': self.openRight,
        })

        return result

    @classmethod
    def fromDict(cls, d):
        """
        Create a new instance from attribute values provided in a dictionary.

        @param d: A C{dict} with keys/values for the attributes of a new
            instance of this class. Keys 'id' and 'sequence' with C{str} values
            must be provided. A 'quality' C{str} key is optional. Keys 'start'
            and 'stop' must have C{int} values. Keys 'openLeft' and 'openRight'
            are C{bool}, all keys are as described in the docstring for this
            class.
        @return: A new instance of this class, with values taken from C{d}.
        """
        # Make a dummy instance whose attributes we can set explicitly.
        new = cls(AARead('', ''), 0, 0, True, True)
        new.id = d['id']
        new.sequence = d['sequence']
        new.quality = d.get('quality')
        new.start = d['start']
        new.stop = d['stop']
        new.openLeft = d['openLeft']
        new.openRight = d['openRight']
        return new


class SSAARead(AARead):
    """
    Hold information to work with AAReads that have secondary structure
    information attached to them.

    Note that this class (currently) has no quality string associated with it.

    @param id: A C{str} describing the read.
    @param sequence: A C{str} of sequence information.
    @param structure: A C{str} of structure information.
    @raise ValueError: If the sequence and structure lengths are not the same.
    """
    def __init__(self, id, sequence, structure):
        if six.PY3:
            super().__init__(id, sequence)
        else:
            AARead.__init__(self, id, sequence)
        self.structure = structure

        if len(sequence) != len(structure):
            raise ValueError(
                'Invalid read: sequence length (%d) != structure length (%d)' %
                (len(sequence), len(structure)))

    def __eq__(self, other):
        return (self.id == other.id and
                self.sequence == other.sequence and
                self.structure == other.structure)

    def __hash__(self):
        """
        Calculate a hash key for a read.

        @return: The C{int} hash key for the read.
        """
        return hash(md5(self.id.encode('UTF-8') + b'\0' +
                        self.sequence.encode('UTF-8') + b'\0' +
                        self.structure.encode('UTF-8')).digest())

    def __getitem__(self, item):
        sequence = self.sequence[item]
        structure = None if self.structure is None else self.structure[item]
        return self.__class__(self.id, sequence, structure)

    def toString(self, format_='fasta-ss', structureSuffix=':structure'):
        """
        Convert the read to a string in PDB format (sequence & structure). This
        consists of two FASTA records, one for the sequence then one for the
        structure.

        @param format_: Either 'fasta-ss' or 'fasta'. In the former case, the
            structure information is returned. Otherwise, plain FASTA is
            returned.
        @param structureSuffix: The C{str} suffix to append to the read id
            for the second FASTA record, containing the structure information.
        @raise ValueError: If C{format_} is not 'fasta'.
        @return: A C{str} representing the read sequence and structure in PDB
            FASTA format.
        """
        if format_ == 'fasta-ss':
            return '>%s\n%s\n>%s%s\n%s\n' % (
                self.id, self.sequence, self.id, structureSuffix,
                self.structure)
        else:
            if six.PY3:
                return super().toString(format_=format_)
            else:
                return AARead.toString(self, format_=format_)

    def toDict(self):
        """
        Get information about this read in a dictionary.

        @return: A C{dict} with keys/values for the attributes of self.
        """
        return {
            'id': self.id,
            'sequence': self.sequence,
            'structure': self.structure,
        }

    @classmethod
    def fromDict(cls, d):
        """
        Create a new instance from attribute values provided in a dictionary.

        @param d: A C{dict} with keys/values for the attributes of a new
            instance of this class. Keys 'id', 'sequence', and 'structure'
            with C{str} values must be provided.
        @return: A new instance of this class, with values taken from C{d}.
        """
        return cls(d['id'], d['sequence'], d['structure'])

    def newFromSites(self, sites, exclude=False):
        """
        Create a new read from self, with only certain sites.

        @param sites: A set of C{int} 0-based sites (i.e., indices) in
            sequences that should be kept. If C{None} (the default), all sites
            are kept.
        @param exclude: If C{True} the C{sites} will be excluded, not
            included.
        """
        if exclude:
            sites = set(range(len(self))) - sites

        newSequence = []
        newStructure = []
        for index, (base, structure) in enumerate(zip(self.sequence,
                                                      self.structure)):
            if index in sites:
                newSequence.append(base)
                newStructure.append(structure)
        read = self.__class__(self.id, ''.join(newSequence),
                              ''.join(newStructure))

        return read


class SSAAReadWithX(SSAARead):
    """
    Hold information and methods to work with C{SSAARead}s allowing 'X'
    characters to appear in sequences.
    """
    ALPHABET = set(AA_LETTERS + ['X'])


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

    def toDict(self):
        """
        Get information about this read in a dictionary.

        @return: A C{dict} with keys/values for the attributes of self.
        """
        if six.PY3:
            result = super().toDict()
        else:
            result = AARead.toDict(self)

        result.update({
            'frame': self.frame,
            'reverseComplemented': self.reverseComplemented,
        })

        return result

    @classmethod
    def fromDict(cls, d):
        """
        Create a new instance from attribute values provided in a dictionary.

        @param d: A C{dict} with keys/values for the attributes of a new
            instance of this class. Keys 'id' and 'sequence' with C{str} values
            must be provided. A 'quality' C{str} key is optional. Key 'frame'
            must have an C{int} value. Key 'reverseComplemented' must be a
            C{bool}, all keys are as described in the docstring for this class.
        @return: A new instance of this class, with values taken from C{d}.
        """
        # Make a dummy instance whose attributes we can set explicitly.
        new = cls(AARead('', ''), 0, True)
        new.id = d['id']
        new.sequence = d['sequence']
        new.quality = d.get('quality')
        new.frame = d['frame']
        new.reverseComplemented = d['reverseComplemented']
        return new

    def maximumORFLength(self, openORFs=True):
        """
        Return the length of the longest (possibly partial) ORF in a translated
        read. The ORF may originate or terminate outside the sequence, which is
        why the length is just a lower bound.
        """
        return max(len(orf) for orf in self.ORFs(openORFs))


class ReadFilter(object):
    """
    Create a function that can be used to filter a set of reads to produce a
    desired subset.

    Note: there are many additional filtering options that could be added,
    e.g., on complexity fraction, on GC %, on quality, etc.

    @param minLength: The minimum acceptable length.
    @param maxLength: The maximum acceptable length.
    @param removeGaps: If C{True} remove all gaps ('-' characters) from the
        read sequences.
    @param whitelist: If not C{None}, a set of exact read ids that are
        always acceptable (though other characteristics, such as length,
        of a whitelisted id may rule it out).
    @param blacklist: If not C{None}, a set of exact read ids that are
        never acceptable.
    @param whitelistFile: If not C{None}, a C{str} filename containing lines
        that give exact ids that are always acceptable.
    @param blacklistFile: If not C{None}, a C{str} filename containing lines
        that give exact ids that are never acceptable.
    @param titleRegex: A regex that read ids must match.
    @param negativeTitleRegex: A regex that read ids must not match.
    @param truncateTitlesAfter: A string that read ids will be truncated
        beyond. If the truncated version of an id has already been seen,
        that sequence will be skipped.
    @param keepSequences: Either C{None} or a set of C{int}s corresponding
        to reads that should be kept. Indexing starts at zero. The numbers
        refer to the sequential number of the sequence in the list of reads.
    @param removeSequences: Either C{None} or a set of C{int}s corresponding
        to reads that should be removed. Indexing starts at zero. The numbers
        refer to the sequential number of the sequence in the list of reads.
    @param head: If not C{None}, the C{int} number of sequences at the
        start of the reads to return. Later sequences are skipped.
    @param removeDuplicates: If C{True} remove duplicated reads based only on
        sequence identity.
    @param removeDuplicatesById: If C{True} remove duplicated reads based
        only on read id.
    @param removeDescriptions: If C{True} remove the description (the part
        following the first whitespace) from read ids. The description is
        removed after applying the function specified by --idLambda (if any).
    @param modifier: If not C{None}, a function that is passed a read
        and which either returns a read or C{None}. If it returns a read,
        that read is passed through the filter. If it returns C{None},
        the read is omitted. Such a function can be used to do customized
        filtering, to change sequence ids, etc.
    @param randomSubset: If not C{None}, an C{int} giving the number of
        sequences that should be returned. These will be selected at
        random, using an algorithm that does a single pass over the data
        and which needs to know (in advance) the total number of reads.
        Because a C{Reads} instance does not always know its length, it is
        possible to specify the number of reads by passing a C{trueLength}
        argument. Note that the random selection is done before any other
        filtering. Due to this, if you want to extract a random subset of
        the reads filtered in another way, it will be best to call filter
        twice rather than doing both types of filtering in one step. E.g.,
        you very likely should do this:
            reads.filter(maxLength=100).filter(randomSubset=20)
        rather than this:
            reads.filter(maxLength=100, randomSubset=20)
        The second version will extract a random subset of 20 reads and
        only return those that have length <= 100, so your result may have
        less than 20 reads. The former version extracts reads of the
        desired length and then takes 20 reads at random from that set, so
        you'll always get 20 reads in your result, assuming there are at
        least that many reads satisfying the length filter.
    @param trueLength: The C{int} number of reads in this C{Reads} instance.
        Under normal circumstances it will not be necessary to pass this
        argument. However in some cases a subclass (e.g., with
        C{dark.fasta.FastaReads}) does not know its length until its data
        has been read from disk. In such cases, it is not possible to
        choose a random subset without keeping the subset in memory (which
        is undesirable). See
        https://en.wikipedia.org/wiki/Reservoir_sampling for one approach
        when the set size is unknown. However, it is possible to filter a
        random subset in a single pass over the data without keeping the
        set in memory if the set size is known. C{trueLength} makes it
        possible to pass the actual number of reads (this will obviously
        need to be obtained via some other mechanism).
    @param sampleFraction: If not C{None}, a [0.0, 1.0] C{float} indicating
        a fraction of the reads that should be allowed to pass through the
        filter. The sample size will only be approximately the product of
        the C{sampleFraction} and the number of reads. The sample is taken
        at random. If you try to combine this filter with C{randomSubset}
        a C{ValueError} will be raised. If you need both filters, run them
        one after another.
    @param sequenceNumbersFile: If not C{None}, gives the C{str} name of a
        file containing (1-based) sequence numbers, in ascending order,
        one per line. Only those sequences matching the given numbers will
        be kept.
    @param idLambda: If not C{None}, a C{str} Python lambda function
        specification to use to modify read ids. The function is applied
        before removing the description (if --removeDescriptions is also
        specified).
    @param readLambda: If not C{None}, a C{str} Python lambda function
        specification to use to modify reads. The function will be passed,
        and must return, a single Read (or one of its subclasses). This
        function is called after the --idLambda function, if any.
    @param keepSites: A set of C{int} 0-based sites (i.e., indices) in
        sequences that should be kept. If C{None} (the default), all sites are
        kept.
    @param removeSites: A set of C{int} 0-based sites (i.e., indices) in
        sequences that should be removed. If C{None} (the default), no sites
        are removed.
    @param reverse: If C{True}, reverse the sequences. Reversing happens
        at a very late stage (i.e., after sites are altered via keepSites
        and removeSites).
    @param reverseComplement: If C{True}, replace seqeunces with their reverse
        complements. Reversing happens at a very late stage (i.e., after sites
        are altered via keepSites and removeSites).
    @raises ValueError: If C{randomSubset} and C{sampleFraction} are both
        specified, or if C{randomSubset} is specified but C{trueLength} is not,
        or if the sequence numbers in C{sequenceNumbersFile} are
        non-positive or not ascending, or if both C{keepSites} and
        C{removeSites} are given, or if both C{keepSequences} and
        C{removeSequences} are given.
    @raise AttributeError: If C{reverseComplement} is C{True} but the read
        type does not allow for reverse complementing.
    """

    # TODO, when/if needed: make it possible to pass a seed for the RNG
    # when randomSubset or sampleFraction are used. Also possible is to
    # save and restore the state of the RNG and/or to optionally add
    # 'seed=XXX' to the end of the id of the first read, etc.

    def __init__(self, minLength=None, maxLength=None, removeGaps=False,
                 whitelist=None, blacklist=None,
                 whitelistFile=None, blacklistFile=None,
                 titleRegex=None, negativeTitleRegex=None,
                 truncateTitlesAfter=None, keepSequences=None,
                 removeSequences=None, head=None,
                 removeDuplicates=False, removeDuplicatesById=False,
                 removeDescriptions=False, modifier=None, randomSubset=None,
                 trueLength=None, sampleFraction=None,
                 sequenceNumbersFile=None, idLambda=None, readLambda=None,
                 keepSites=None, removeSites=None, reverse=False,
                 reverseComplement=False):

        if randomSubset is not None:
            if sampleFraction is not None:
                raise ValueError('randomSubset and sampleFraction cannot be '
                                 'used simultaneously in a filter. Make two '
                                 'read filters instead.')

            if trueLength is None:
                raise ValueError('trueLength must be supplied if randomSubset '
                                 'is specified.')

        self.minLength = minLength
        self.maxLength = maxLength
        self.removeGaps = removeGaps
        self.head = head
        self.removeDuplicates = removeDuplicates
        self.removeDuplicatesById = removeDuplicatesById
        self.removeDescriptions = removeDescriptions
        self.modifier = modifier
        self.randomSubset = randomSubset
        self.trueLength = trueLength

        if keepSequences and removeSequences:
            raise ValueError(
                'Cannot simultaneously filter using keepSequences and '
                'removeSequences. Call filter twice in succession instead.')
        self.keepSequences = keepSequences
        self.removeSequences = removeSequences

        if keepSites and removeSites:
            raise ValueError(
                'Cannot simultaneously filter using keepSites and '
                'removeSites. Call filter twice in succession instead.')
        self.keepSites = keepSites
        self.removeSites = removeSites

        if reverseComplement:
            # Make sure reverse is not also set.
            reverse = False
        self.reverse = reverse
        self.reverseComplement = reverseComplement

        self.alwaysFalse = False
        self.yieldCount = 0
        self.readIndex = -1

        def _wantedSequences(filename):
            """
            Read and yield integer sequence numbers from a file.

            @raise ValueError: If the sequence numbers are not all positive or
                are not ascending.
            @return: A generator that yields C{int} sequence numbers.
            """
            with open(filename) as fp:
                lastNumber = None
                for line in fp:
                    n = int(line)
                    if lastNumber is None:
                        if n < 1:
                            raise ValueError(
                                'First line of sequence number file %r must '
                                'be at least 1.' % filename)
                        lastNumber = n
                        yield n
                    else:
                        if n > lastNumber:
                            lastNumber = n
                            yield n
                        else:
                            raise ValueError(
                                'Line number file %r contains non-ascending '
                                'numbers %d and %d.' %
                                (filename, lastNumber, n))

        self.wantedSequenceNumberGeneratorExhausted = False
        self.nextWantedSequenceNumber = None

        if sequenceNumbersFile is not None:
            self.wantedSequenceNumberGenerator = _wantedSequences(
                sequenceNumbersFile)
            try:
                self.nextWantedSequenceNumber = next(
                    self.wantedSequenceNumberGenerator)
            except StopIteration:
                # There was a sequence number file, but it was empty. So no
                # reads will ever be accepted.
                self.alwaysFalse = True

        if (whitelist or blacklist or whitelistFile or blacklistFile or
                titleRegex or negativeTitleRegex or truncateTitlesAfter):
            self.titleFilter = TitleFilter(
                whitelist=whitelist, blacklist=blacklist,
                whitelistFile=whitelistFile, blacklistFile=blacklistFile,
                positiveRegex=titleRegex, negativeRegex=negativeTitleRegex,
                truncateAfter=truncateTitlesAfter)
        else:
            self.titleFilter = None

        if removeDuplicates:
            self.sequencesSeen = set()

        if removeDuplicatesById:
            self.idsSeen = set()

        if sampleFraction is not None:
            if sampleFraction == 0.0:
                # The filter method should always return False.
                self.alwaysFalse = True
            elif sampleFraction == 1.0:
                # Passing 1.0 can be treated the same as passing no value.
                # This makes the filter code below simpler.
                sampleFraction = None
        self.sampleFraction = sampleFraction

        self.idLambda = eval(idLambda) if idLambda else None
        self.readLambda = eval(readLambda) if readLambda else None

    def filter(self, read):
        """
        Check if a read passes the filter.

        @param read: A C{Read} instance.
        @return: C{read} if C{read} passes the filter, C{False} if not.
        """
        self.readIndex += 1

        if self.alwaysFalse:
            return False

        if self.wantedSequenceNumberGeneratorExhausted:
            return False

        if self.nextWantedSequenceNumber is not None:
            if self.readIndex + 1 == self.nextWantedSequenceNumber:
                # We want this sequence.
                try:
                    self.nextWantedSequenceNumber = next(
                        self.wantedSequenceNumberGenerator)
                except StopIteration:
                    # The sequence number iterator ran out of sequence
                    # numbers.  We must let the rest of the filtering
                    # continue for the current sequence in case we
                    # throw it out for other reasons (as we might have
                    # done for any of the earlier wanted sequence
                    # numbers).
                    self.wantedSequenceNumberGeneratorExhausted = True
            else:
                # This sequence isn't one of the ones that's wanted.
                return False

        if (self.sampleFraction is not None and
                uniform(0.0, 1.0) > self.sampleFraction):
            # Note that we don't have to worry about the 0.0 or 1.0
            # cases in the above 'if', as they have been dealt with
            # in self.__init__.
            return False

        if self.randomSubset is not None:
            if self.yieldCount == self.randomSubset:
                # The random subset has already been fully returned.
                # There's no point in going any further through the input.
                self.alwaysFalse = True
                return False
            elif uniform(0.0, 1.0) > ((self.randomSubset - self.yieldCount) /
                                      (self.trueLength - self.readIndex)):
                return False

        if self.head is not None and self.readIndex == self.head:
            # We're completely done.
            self.alwaysFalse = True
            return False

        readLen = len(read)
        if ((self.minLength is not None and readLen < self.minLength) or
                (self.maxLength is not None and readLen > self.maxLength)):
            return False

        if self.removeGaps:
            if read.quality is None:
                read = read.__class__(read.id, read.sequence.replace('-', ''))
            else:
                newSequence = []
                newQuality = []
                for base, quality in zip(read.sequence, read.quality):
                    if base != '-':
                        newSequence.append(base)
                        newQuality.append(quality)
                read = read.__class__(
                    read.id, ''.join(newSequence), ''.join(newQuality))

        if (self.titleFilter and
                self.titleFilter.accept(read.id) == TitleFilter.REJECT):
            return False

        if (self.keepSequences is not None and
                self.readIndex not in self.keepSequences):
            return False

        if (self.removeSequences is not None and
                self.readIndex in self.removeSequences):
            return False

        if self.removeDuplicates:
            if read.sequence in self.sequencesSeen:
                return False
            self.sequencesSeen.add(read.sequence)

        if self.removeDuplicatesById:
            if read.id in self.idsSeen:
                return False
            self.idsSeen.add(read.id)

        if self.modifier:
            modified = self.modifier(read)
            if modified is None:
                return False
            else:
                read = modified

        # We have to use 'is not None' in the following tests so the empty set
        # is processed properly.
        if self.keepSites is not None:
            read = read.newFromSites(self.keepSites)
        elif self.removeSites is not None:
            read = read.newFromSites(self.removeSites, exclude=True)

        if self.idLambda:
            newId = self.idLambda(read.id)
            if newId is None:
                return False
            else:
                read.id = newId

        if self.readLambda:
            newRead = self.readLambda(read)
            if newRead is None:
                return False
            else:
                read = newRead

        if self.removeDescriptions:
            read.id = read.id.split()[0]

        if self.reverse:
            read = read.reverse()
        elif self.reverseComplement:
            read = read.reverseComplement()

        self.yieldCount += 1
        return read


# Provide a mapping from all read class names to read classes. This can be
# useful in deserialization.
readClassNameToClass = {
    'AARead': AARead,
    'AAReadORF': AAReadORF,
    'AAReadWithX': AAReadWithX,
    'DNARead': DNARead,
    'RNARead': RNARead,
    'Read': Read,
    'SSAARead': SSAARead,
    'SSAAReadWithX': SSAAReadWithX,
    'TranslatedRead': TranslatedRead,
}

_DNA = set('ACGT')
_RNA = set('ACGU')
_AA = set(AA_NAMES)

# Map read class names to the set of unambiguous sequence bases (loosely
# speaking).
unambiguousBases = {
    'AARead': _AA,
    'AAReadORF': _AA,
    'AAReadWithX': _AA,
    'DNARead': _DNA,
    'RNARead': _RNA,
    'Read': _DNA,
    'SSAARead': _AA,
    'SSAAReadWithX': _AA,
    'TranslatedRead': _AA,
}


class Reads(object):
    """
    Maintain a collection of sequence reads.

    @param initialReads: If not C{None}, an iterable of C{Read} (or C{Read}
        subclass) instances.
    """

    def __init__(self, initialReads=None):
        self._initialReads = initialReads
        self._additionalReads = []
        self._filters = []
        self._iterated = False

    def filterRead(self, read):
        """
        Filter a read, according to our set of filters.

        @param read: A C{Read} instance or one of its subclasses.
        @return: C{False} if the read fails any of our filters, else the
            C{Read} instance returned by our list of filters.
        """
        for filterFunc in self._filters:
            filteredRead = filterFunc(read)
            if filteredRead is False:
                return False
            else:
                read = filteredRead
        return read

    def add(self, read):
        """
        Add a read to this collection of reads.

        @param read: A C{Read} instance.
        """
        self._additionalReads.append(read)

    def __iter__(self):
        """
        Iterate through all the reads.

        @return: A generator that yields reads. The returned read types depend
            on the kind of reads that were added to this instance.
        """
        # self._additionalReads is a regular list.
        for read in self._additionalReads:
            filteredRead = self.filterRead(read)
            if filteredRead is not False:
                yield filteredRead

        _unfilteredLength = len(self._additionalReads)

        # self._initialReads may be a Reads instance and/or may not support
        # len().
        initialReads = self._initialReads or []
        initialReadsLength = 0
        for read in initialReads:
            initialReadsLength += 1
            filteredRead = self.filterRead(read)
            if filteredRead is not False:
                yield filteredRead

        if isinstance(initialReads, Reads):
            _unfilteredLength += initialReads.unfilteredLength()
        else:
            _unfilteredLength += initialReadsLength

        # The value returned by self.iter() may be a Reads instance and/or
        # may not support len().
        subclassReads = self.iter()
        subclassReadsLength = 0
        for read in subclassReads:
            subclassReadsLength += 1
            filteredRead = self.filterRead(read)
            if filteredRead is not False:
                yield filteredRead

        if isinstance(subclassReads, Reads):
            _unfilteredLength += subclassReads.unfilteredLength()
        else:
            _unfilteredLength += subclassReadsLength

        self._unfilteredLength = _unfilteredLength
        self._iterated = True

    def unfilteredLength(self):
        """
        Return the underlying number of reads in C{self}, irrespective of any
        filtering that has been applied.

        To obtain the number of reads in a filtered C{Reads} instance, you
        must count the reads yourself as you iterate it.

        @raises RuntimeError: If C{self} has not been fully iterated.
        @return: The C{int} number of reads in C{self}.
        """
        if self._iterated:
            return self._unfilteredLength
        else:
            raise RuntimeError(
                'The unfiltered length of a Reads instance is unknown until '
                'it has been iterated.')

    def iter(self):
        """
        Placeholder to allow subclasses to provide reads.

        These might be extracted from a file. E.g., the
        C{dark.reads.fasta.FastaReads} class (a subclass of C{Reads})
        overrides this method to provide reads from a file.

        @return: An iterable of C{Read} instances.
        """
        return []

    def save(self, filename, format_='fasta'):
        """
        Write the reads to C{filename} in the requested format.

        @param filename: Either a C{str} file name to save into (the file will
            be overwritten) or an open file descriptor (e.g., sys.stdout).
        @param format_: A C{str} format to save as, either 'fasta', 'fastq' or
            'fasta-ss'.
        @raise ValueError: if C{format_} is 'fastq' and a read with no quality
            is present, or if an unknown format is requested.
        @return: An C{int} giving the number of reads in C{self}.
        """
        format_ = format_.lower()
        count = 0

        if isinstance(filename, str):
            try:
                with open(filename, 'w') as fp:
                    for read in self:
                        fp.write(read.toString(format_))
                        count += 1
            except ValueError:
                unlink(filename)
                raise
        else:
            # We have a file-like object.
            for read in self:
                filename.write(read.toString(format_))
                count += 1
        return count

    def filter(self, **kwargs):
        """
        Add a filter to this C{Reads} instance.

        @param kwargs: Keyword arguments, as accepted by C{ReadFilter}.
        @return: C{self}.
        """
        readFilter = ReadFilter(**kwargs)
        self._filters.append(readFilter.filter)
        return self

    def clearFilters(self):
        """
        Clear all filters on this C{Reads} instance.

        @return: C{self}.
        """
        self._filters = []
        return self

    def summarizePosition(self, index):
        """
        Compute residue counts at a specific sequence index.

        @param index: an C{int} index into the sequence.
        @return: A C{dict} with the count of too-short (excluded) sequences,
            and a Counter instance giving the residue counts.
        """
        countAtPosition = Counter()
        excludedCount = 0

        for read in self:
            try:
                countAtPosition[read.sequence[index]] += 1
            except IndexError:
                excludedCount += 1

        return {
            'excludedCount': excludedCount,
            'countAtPosition': countAtPosition
        }

    def sitesMatching(self, targets, matchCase, any_):
        """
        Find sites (i.e., sequence indices) that match a given set of target
        sequence bases.

        @param targets: A C{set} of sequence bases to look for.
        @param matchCase: If C{True}, case will be considered in matching.
        @param any_: If C{True}, return sites that match in any read. Else
            return sites that match in all reads.
        @return: A C{set} of 0-based sites that indicate where the target
            bases occur in our reads. An index will be in this set if any of
            our reads has any of the target bases in that location.
        """
        # If case is unimportant, we convert everything (target bases and
        # sequences, as we read them) to lower case.
        if not matchCase:
            targets = set(map(str.lower, targets))

        result = set() if any_ else None
        for read in self:
            sequence = read.sequence if matchCase else read.sequence.lower()
            matches = set(index for (index, base) in enumerate(sequence)
                          if base in targets)
            if any_:
                result |= matches
            else:
                if result is None:
                    result = matches
                else:
                    result &= matches
                # We can exit early if we run out of possible sites.
                if not result:
                    break

        # Make sure we don't return None.
        return result or set()


class ReadsInRAM(Reads):
    """
    Maintain a collection of sequence reads in RAM.

    @param initialReads: If not C{None}, an iterable of C{Read} (or a C{Read}
        subclass) instances.
    """

    # This class provides some C{list} like methods (len and indexing) but
    # is not an actual list or list subclass. That's because we want to inherit
    # the methods of C{Reads}, and I considered it too messy to use double
    # inheritance. If you want a real list, you can just call C{list} on a
    # C{Reads} or C{ReadsInRAM} instance.

    def __init__(self, initialReads=None):
        if six.PY3:
            super().__init__(initialReads)
        else:
            Reads.__init__(self, initialReads)

        # Read all initial reads into memory.
        if initialReads:
            for read in initialReads:
                self.add(read)

        # Set self._iterated to True in case someone calls unfilteredLength
        # (see Reads).
        self._iterated = True

    def __len__(self):
        return self._additionalReads.__len__()

    def __getitem__(self, item):
        return self._additionalReads.__getitem__(item)

    def __setitem__(self, item, value):
        return self._additionalReads.__setitem__(item, value)

    def __iter__(self):
        return self._additionalReads.__iter__()


def addFASTACommandLineOptions(parser):
    """
    Add standard command-line options to an argparse parser.

    @param parser: An C{argparse.ArgumentParser} instance.
    """

    parser.add_argument(
        '--fastaFile', type=open, default=sys.stdin, metavar='FILENAME',
        help=('The name of the FASTA input file. Standard input will be read '
              'if no file name is given.'))

    parser.add_argument(
        '--readClass', default='DNARead', choices=readClassNameToClass,
        metavar='CLASSNAME',
        help=('If specified, give the type of the reads in the input. '
              'Possible choices: %s.' % ', '.join(readClassNameToClass)))

    # A mutually exclusive group for either --fasta, --fastq, or --fasta-ss
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        '--fasta', default=False, action='store_true',
        help=('If specified, input will be treated as FASTA. This is the '
              'default.'))

    group.add_argument(
        '--fastq', default=False, action='store_true',
        help='If specified, input will be treated as FASTQ.')

    group.add_argument(
        '--fasta-ss', dest='fasta_ss', default=False, action='store_true',
        help=('If specified, input will be treated as PDB FASTA '
              '(i.e., regular FASTA with each sequence followed by its '
              'structure).'))


def parseFASTACommandLineOptions(args):
    """
    Examine parsed command-line options and return a Reads instance.

    @param args: An argparse namespace, as returned by the argparse
        C{parse_args} function.
    @return: A C{Reads} subclass instance, depending on the type of FASTA file
        given.
    """
    # Set default FASTA type.
    if not (args.fasta or args.fastq or args.fasta_ss):
        args.fasta = True

    readClass = readClassNameToClass[args.readClass]

    if args.fasta:
        from dark.fasta import FastaReads
        return FastaReads(args.fastaFile, readClass=readClass)
    elif args.fastq:
        from dark.fastq import FastqReads
        return FastqReads(args.fastaFile, readClass=readClass)
    else:
        from dark.fasta_ss import SSFastaReads
        return SSFastaReads(args.fastaFile, readClass=readClass)

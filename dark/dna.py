from __future__ import division

from collections import defaultdict

from dark.utils import countPrint

try:
    from itertools import zip_longest
except ImportError:
    # zip_longest does not exist in Python 2.7 itertools. We should be able
    # to get it via from six.moves import zip_longest according to
    # https://pythonhosted.org/six/index.html?highlight=zip_longest but
    # that doesn't work for me.
    from itertools import izip_longest as zip_longest

# A list of the ambiguous values is given at
# https://en.wikipedia.org/wiki/Nucleic_acid_notation
AMBIGUOUS = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'M': {'A', 'C'},
    'R': {'A', 'G'},
    'W': {'A', 'T'},
    'S': {'G', 'C'},
    'K': {'G', 'T'},
    'Y': {'C', 'T'},
    'V': {'A', 'C', 'G'},
    'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'},
    'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'},
}

# Make a reverse version of AMBIGUOUS.
BASES_TO_AMBIGUOUS = dict(
    (''.join(sorted(bases)), symbol) for symbol, bases in AMBIGUOUS.items())


def matchToString(dnaMatch, read1, read2, matchAmbiguous=True, indent='',
                  offsets=None, includeGapLocations=True):
    """
    Format a DNA match as a string.

    @param dnaMatch: A C{dict} returned by C{compareDNAReads}.
    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param matchAmbiguous: If C{True}, ambiguous nucleotides that are
        possibly correct were counted as actually being correct. Otherwise,
        the match was done strictly, insisting that only non-ambiguous
        nucleotides could contribute to the matching nucleotide count.
    @param indent: A C{str} to indent all returned lines with.
    @param offsets: If not C{None}, a C{set} of offsets of interest that were
        only considered when making C{match}.
    @param includeGapLocations: If C{True} indicate the (1-based) locations of
        gaps.
    @return: A C{str} describing the match.
    """
    match = dnaMatch['match']
    identicalMatchCount = match['identicalMatchCount']
    ambiguousMatchCount = match['ambiguousMatchCount']
    gapMismatchCount = match['gapMismatchCount']
    gapGapMismatchCount = match['gapGapMismatchCount']
    nonGapMismatchCount = match['nonGapMismatchCount']

    if offsets:
        len1 = len2 = len(offsets)
    else:
        len1, len2 = map(len, (read1, read2))

    result = []
    append = result.append

    append(countPrint('%sExact matches' % indent, identicalMatchCount,
                      len1, len2))
    append(countPrint('%sAmbiguous matches' % indent, ambiguousMatchCount,
                      len1, len2))
    if ambiguousMatchCount and identicalMatchCount:
        anyMatchCount = identicalMatchCount + ambiguousMatchCount
        append(countPrint('%sExact or ambiguous matches' % indent,
                          anyMatchCount, len1, len2))

    mismatchCount = (gapMismatchCount + gapGapMismatchCount +
                     nonGapMismatchCount)
    append(countPrint('%sMismatches' % indent, mismatchCount, len1, len2))
    conflicts = 'conflicts' if matchAmbiguous else 'conflicts or ambiguities'
    append(countPrint('%s  Not involving gaps (i.e., %s)' % (indent,
                      conflicts), nonGapMismatchCount, len1, len2))
    append(countPrint('%s  Involving a gap in one sequence' % indent,
                      gapMismatchCount, len1, len2))
    append(countPrint('%s  Involving a gap in both sequences' % indent,
                      gapGapMismatchCount, len1, len2))

    for read, key in zip((read1, read2), ('read1', 'read2')):
        append('%s  Id: %s' % (indent, read.id))
        length = len(read)
        append('%s    Length: %d' % (indent, length))
        gapCount = len(dnaMatch[key]['gapOffsets'])
        append(countPrint('%s    Gaps' % indent, gapCount, length))
        if includeGapLocations and gapCount:
            append(
                '%s    Gap locations (1-based): %s' %
                (indent,
                 ', '.join(map(lambda offset: str(offset + 1),
                               sorted(dnaMatch[key]['gapOffsets'])))))
        ambiguousCount = len(dnaMatch[key]['ambiguousOffsets'])
        append(countPrint('%s    Ambiguous' % indent, ambiguousCount, length))
        extraCount = dnaMatch[key]['extraCount']
        if extraCount:
            append(countPrint('%s    Extra nucleotides at end' % indent,
                              extraCount, length))

    return '\n'.join(result)


def compareDNAReads(read1, read2, matchAmbiguous=True, gapChars='-',
                    offsets=None):
    """
    Compare two DNA sequences.

    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct, and score these in the
        ambiguousMatchCount. Otherwise, we are strict and insist that only
        non-ambiguous nucleotides can contribute to the matching nucleotide
        count.
    @param gapChars: An object supporting __contains__ with characters that
        should be considered to be gaps.
    @param offsets: If not C{None}, a C{set} of offsets of interest. Offsets
        not in the set will not be considered.
    @return: A C{dict} with information about the match and the individual
        sequences (see below).
    """
    identicalMatchCount = ambiguousMatchCount = 0
    gapMismatchCount = nonGapMismatchCount = gapGapMismatchCount = 0
    read1ExtraCount = read2ExtraCount = 0
    read1GapOffsets = []
    read2GapOffsets = []
    read1AmbiguousOffsets = []
    read2AmbiguousOffsets = []
    empty = set()

    def _identicalMatch(a, b):
        return a == b and len(AMBIGUOUS[a]) == 1

    def _ambiguousMatch(a, b, matchAmbiguous):
        """
        Checks if two characters match ambiguously if matchAmbiguous is True.
        A match is an ambiguous match if it is not an identical match, but the
        sets of ambiguous characters overlap.
        """
        return (matchAmbiguous and
                not _identicalMatch(a, b) and
                AMBIGUOUS.get(a, empty) & AMBIGUOUS.get(b, empty))

    for offset, (a, b) in enumerate(zip_longest(read1.sequence.upper(),
                                                read2.sequence.upper())):
        # Use 'is not None' in the following to allow an empty offsets set
        # to be passed.
        if offsets is not None and offset not in offsets:
            continue
        if len(AMBIGUOUS.get(a, '')) > 1:
            read1AmbiguousOffsets.append(offset)
        if len(AMBIGUOUS.get(b, '')) > 1:
            read2AmbiguousOffsets.append(offset)
        if a is None:
            # b has an extra character at its end (it cannot be None).
            assert b is not None
            read2ExtraCount += 1
            if b in gapChars:
                read2GapOffsets.append(offset)
        elif b is None:
            # a has an extra character at its end.
            read1ExtraCount += 1
            if a in gapChars:
                read1GapOffsets.append(offset)
        else:
            # We have a character from both sequences (they could still be
            # gap characters).
            if a in gapChars:
                read1GapOffsets.append(offset)
                if b in gapChars:
                    # Both are gaps. This can happen (though hopefully not
                    # if the sequences were pairwise aligned).
                    gapGapMismatchCount += 1
                    read2GapOffsets.append(offset)
                else:
                    # a is a gap, b is not.
                    gapMismatchCount += 1
            else:
                if b in gapChars:
                    # b is a gap, a is not.
                    gapMismatchCount += 1
                    read2GapOffsets.append(offset)
                else:
                    # Neither is a gap character.
                    if _identicalMatch(a, b):
                        identicalMatchCount += 1
                    elif _ambiguousMatch(a, b, matchAmbiguous):
                        ambiguousMatchCount += 1
                    else:
                        nonGapMismatchCount += 1

    return {
        'match': {
            'identicalMatchCount': identicalMatchCount,
            'ambiguousMatchCount': ambiguousMatchCount,
            'gapMismatchCount': gapMismatchCount,
            'gapGapMismatchCount': gapGapMismatchCount,
            'nonGapMismatchCount': nonGapMismatchCount,
        },
        'read1': {
            'ambiguousOffsets': read1AmbiguousOffsets,
            'extraCount': read1ExtraCount,
            'gapOffsets': read1GapOffsets,
        },
        'read2': {
            'ambiguousOffsets': read2AmbiguousOffsets,
            'extraCount': read2ExtraCount,
            'gapOffsets': read2GapOffsets,
        },
    }


def findKozakConsensus(read):
    """
    In a given DNA sequence, search for a Kozak consensus: (gcc)gccRccATGG.
    The upper case bases in that pattern are required, and the lower case
    bases are the ones most frequently found at the given positions. The
    initial 'gcc' sequence (in parentheses) is of uncertain significance
    and is not taken into account here.

    @param read: A C{DNARead} instance to be checked for Kozak consensi.
    @return: A generator that yields C{DNAKozakRead} instances.
    """
    from dark.reads import DNAKozakRead

    readLen = len(read)
    if readLen > 9:
        offset = 6
        readSeq = read.sequence
        while offset < readLen - 3:
            triplet = readSeq[offset:offset + 3]
            if triplet == 'ATG':
                if readSeq[offset + 3] == 'G':
                    if readSeq[offset - 3] in 'GA':
                        kozakQualityCount = sum((
                            readSeq[offset - 1] == 'C',
                            readSeq[offset - 2] == 'C',
                            readSeq[offset - 4] == 'C',
                            readSeq[offset - 5] == 'C',
                            readSeq[offset - 6] == 'G'))

                        kozakQualityPercent = kozakQualityCount / 5.0 * 100
                        yield DNAKozakRead(read, offset - 6, offset + 4,
                                           kozakQualityPercent)
            offset += 1


class FloatBaseCounts(object):
    """
    Hold a floating point count of possible nucleotide bases.

    @param codes: An iterable of nucleotide codes.
    @param unknownAreAmbiguous: If C{True}, any unknown character (e.g., a '-'
        gap or '?' unknown base) will be treated as being fully ambiguous
        (i.e., could be any of ACGT). Otherwise, all unknown characters are
        collected under the count for '-'.
    """
    def __init__(self, codes, unknownAreAmbiguous=False):
        self.codes = list(map(str.upper, codes))
        self.unknownAreAmbiguous = unknownAreAmbiguous
        self.n = len(self.codes)
        counts = defaultdict(float)

        default = self._default = set('ACGT') if unknownAreAmbiguous else {'-'}

        for code in self.codes:
            possible = AMBIGUOUS.get(code, default)
            frac = 1.0 / len(possible)
            for p in possible:
                counts[p] += frac

        # Sort first by base.
        _sorted = [(base, counts[base]) for base in sorted(counts)]

        # Then by count (reversed).
        def key(item):
            return item[1]

        self._sorted = sorted(_sorted, key=key, reverse=True)
        self.counts = counts

    def mostFrequent(self):
        """
        Which bases are the most frequent?

        @return: A C{set} of the most frequent bases.
        """
        maxCount = self._sorted[0][1]
        return set(base for base, count in self._sorted if count == maxCount)

    def highestFrequency(self):
        """
        What is the frequency of the most frequent base?

        @return: The C{float} frequency of the most common base.
        """
        if len(self.counts) < 2:
            return 1.0
        else:
            return self._sorted[0][1] / float(self.n)

    def homogeneous(self, level):
        """
        Does the most frequent base occurs at least C{level} fraction of the
        time?

        @param level: A C{float} fraction.
        @return: C{True} if the most common base occurs at least C{level}
            fraction of the time. If there are no bases at all, this is
            considered homogeneous.
        """
        return self.highestFrequency() >= level

    def __len__(self):
        return len(self.counts)

    def __str__(self):
        fmt = '%d' if all(c == int(c) for b, c in self._sorted) else '%.2f'
        return '%s (%.3f)' % (
            ' '.join(('%s:' + fmt) % (b, c) for b, c in self._sorted),
            self.highestFrequency())

    def variable(self, confirm=True):
        """
        Are the nucleotides variable?

        @param confirm: If C{True}, confirm that there is variability by
            looking at the ambiguous nucleotides. Else just report C{True}
            if there is more than one code (which might not indicate actual
            variability, since two codes could be ambiguous and have a
            nucleotide in common).
        """
        codes = self.codes

        if confirm:
            unambiguous = set()
            ambiguousIntersection = None
            default = self._default
            for code in codes:
                possible = AMBIGUOUS.get(code, default)
                if len(possible) == 1:
                    unambiguous.add(code)
                else:
                    if ambiguousIntersection is None:
                        ambiguousIntersection = set(possible)
                    else:
                        ambiguousIntersection.intersection_update(possible)

            if len(unambiguous) == 0:
                # There were no unambiguous nucleotide codes.

                # Sanity check: there must have been some ambiguous sites.
                assert ambiguousIntersection is not None

                if len(ambiguousIntersection) == 0:
                    # The ambiguous sites had nothing in common, so
                    # variation must exist (it cannot be determined what
                    # the variation is, but we don't care about that).
                    return True
                else:
                    # All the ambiguous sites have at least one nucleotide
                    # in common, so we can't be sure there's any variation.
                    pass
            elif len(unambiguous) == 1:
                # All the unambiguous sites agree. Do any of the ambiguous
                # sites (if any) not allow the unambiguous nucleotide in
                # their possibilities?
                if ambiguousIntersection is None:
                    # There were no ambiguous sites, so there's no
                    # variation here.
                    pass
                else:
                    # If any of the ambiguous sites excludes the single
                    # unambiguous nucleotide, then variation must exist.
                    nt = unambiguous.pop()
                    for code in codes:
                        possible = AMBIGUOUS.get(code, default)
                        if nt not in possible:
                            return True
            elif len(unambiguous) > 1:
                return True
        else:
            if len(codes) > 1:
                return True

        return False


def sequenceToRegex(sequence, wildcards='-?'):
    """
    Convert a potentially ambiguous DNA sequence into a regular expression.
    '?' and '-' are translated into [ACGT].

    @param sequence: A C{str} DNA sequence, possibly with ambiguous codes.
        Case insensitive.
    @param wildcards: A C{set} (or C{str}) with characters that should be
        translated to '[ACGT]'. Note that this happens only if the standard
        ambiguous lookup fails (the order could be changed one day if we need
        to override, or we could allow the passing of an ambiguity mapping).
        Wildcards are case sensitive.
    @raise KeyError: If any character in C{sequence} is unknown.
    @return: A C{str} regular expression with [...] for the ambiguous codes in
        C{sequence}.
    """
    result = []
    append = result.append
    for base in sequence.upper():
        try:
            possible = ''.join(sorted(AMBIGUOUS[base]))
        except KeyError:
            if base in wildcards:
                possible = 'ACGT'
            else:
                raise

        append(('[%s]' % possible) if len(possible) > 1 else possible)

    return ''.join(result)

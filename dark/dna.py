from __future__ import division

from dark.utils import countPrint
from dark.reads import DNAKozakRead
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
                  offsets=None):
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
        if gapCount:
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

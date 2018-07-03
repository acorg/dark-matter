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


def compareDNAReads(read1, read2, matchAmbiguous=True, gapChars=('-'),
                    offsets=None):
    """
    Compare two DNA sequences.

    Note that including 'N' in the default gapChars results (obviously) in
    treating N as a gap. Although you could argue that those two things are
    not the same, for all practical purposes (or at least in the code below)
    they can be treated as if they are. We treat a '?' as a gap character
    too, as it is also completely ambiguous.  In the returned result, the
    gap information (mismatches, indices, etc) therefore include places where
    the sequences have an 'N'. If you don't want 'N' treated in this way, just
    pass a C{gapChars} value without it.

    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct, and score these in the
        ambiguousMatchCount. Otherwise, we are strict and insist that only
        non-ambiguous nucleotides can contribute to the matching nucleotide
        count.
    @param gapChars: An iterable containing characters that should be
        considered to be gaps.
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

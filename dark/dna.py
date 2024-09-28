from collections import defaultdict
from operator import itemgetter
from itertools import zip_longest

from dark.utils import countPrint

# A list of the ambiguous values is given at
# https://en.wikipedia.org/wiki/Nucleic_acid_notation
AMBIGUOUS: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "M": {"A", "C"},
    "R": {"A", "G"},
    "W": {"A", "T"},
    "S": {"G", "C"},
    "K": {"G", "T"},
    "Y": {"C", "T"},
    "V": {"A", "C", "G"},
    "H": {"A", "C", "T"},
    "D": {"A", "G", "T"},
    "B": {"C", "G", "T"},
    "N": {"A", "C", "G", "T"},
}

# Make a reverse version of AMBIGUOUS.
BASES_TO_AMBIGUOUS = dict(
    ("".join(sorted(bases)), symbol) for symbol, bases in AMBIGUOUS.items()
)


def matchToString(
    dnaMatch,
    read1,
    read2,
    matchAmbiguous=True,
    indent="",
    offsets=None,
    includeGapLocations=True,
    includeNoCoverageLocations=True,
    includeAmbiguousMatches=False,
    includeNonGapMismatches=False,
):
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
    @param includeNoCoverageLocations: If C{True} indicate the (1-based)
        locations of no coverage.
    @param includeAmbiguousMatches: If C{True} indicate the (1-based)
        locations of ambiguous matches.
    @param includeNonGapMismatches: If C{True} indicate the (1-based) locations
        of non-gap mismatches.
    @return: A C{str} describing the match.
    """
    match = dnaMatch["match"]
    identicalMatchCount = match["identicalMatchCount"]
    ambiguousMatchCount = match["ambiguousMatchCount"]
    gapMismatchCount = match["gapMismatchCount"]
    gapGapMismatchCount = match["gapGapMismatchCount"]
    nonGapMismatchCount = match["nonGapMismatchCount"]
    noCoverageCount = match["noCoverageCount"]
    noCoverageNoCoverageCount = match["noCoverageNoCoverageCount"]

    if offsets:
        len1 = len2 = len(offsets)
    else:
        len1, len2 = len(read1), len(read2)

    result = []
    append = result.append

    append(countPrint("%sExact matches" % indent, identicalMatchCount, len1, len2))
    append(countPrint("%sAmbiguous matches" % indent, ambiguousMatchCount, len1, len2))

    if noCoverageCount:
        append(
            countPrint(
                "%sExact matches (ignoring no coverage sites)" % indent,
                identicalMatchCount,
                len1 - noCoverageCount,
                len2 - noCoverageCount,
            )
        )
        append(
            countPrint(
                "%sAmbiguous matches (ignoring no coverage sites)" % indent,
                ambiguousMatchCount,
                len1 - noCoverageCount,
                len2 - noCoverageCount,
            )
        )

    if ambiguousMatchCount and identicalMatchCount:
        anyMatchCount = identicalMatchCount + ambiguousMatchCount
        append(
            countPrint(
                "%sExact or ambiguous matches" % indent, anyMatchCount, len1, len2
            )
        )
        if noCoverageCount:
            append(
                countPrint(
                    "%sExact or ambiguous matches (ignoring no coverage sites)"
                    % indent,
                    anyMatchCount,
                    len1 - noCoverageCount,
                    len2 - noCoverageCount,
                )
            )

    mismatchCount = gapMismatchCount + gapGapMismatchCount + nonGapMismatchCount
    append(countPrint("%sMismatches" % indent, mismatchCount, len1, len2))
    conflicts = "conflicts" if matchAmbiguous else "conflicts or ambiguities"
    append(
        countPrint(
            "%s  Not involving gaps (i.e., %s)" % (indent, conflicts),
            nonGapMismatchCount,
            len1,
            len2,
        )
    )
    append(
        countPrint(
            "%s  Involving a gap in one sequence" % indent, gapMismatchCount, len1, len2
        )
    )
    append(
        countPrint(
            "%s  Involving a gap in both sequences" % indent,
            gapGapMismatchCount,
            len1,
            len2,
        )
    )
    append(
        countPrint(
            "%s  Involving no coverage in one sequence" % indent,
            noCoverageCount,
            len1,
            len2,
        )
    )
    append(
        countPrint(
            "%s  Involving no coverage in both sequences" % indent,
            noCoverageNoCoverageCount,
            len1,
            len2,
        )
    )

    for read, key in zip((read1, read2), ("read1", "read2")):
        append("%s  Id: %s" % (indent, read.id))
        length = len(read)
        append("%s    Length: %d" % (indent, length))
        gapCount = len(dnaMatch[key]["gapOffsets"])
        append(countPrint("%s    Gaps" % indent, gapCount, length))
        if includeGapLocations and gapCount:
            append(
                "%s    Gap locations (1-based): %s"
                % (
                    indent,
                    ", ".join(
                        map(
                            lambda offset: str(offset + 1),
                            sorted(dnaMatch[key]["gapOffsets"]),
                        )
                    ),
                )
            )
        ncCount = len(dnaMatch[key]["noCoverageOffsets"])
        append(countPrint("%s    No coverage" % indent, ncCount, length))
        if includeNoCoverageLocations and ncCount:
            append(
                "%s    No coverage locations (1-based): %s"
                % (
                    indent,
                    ", ".join(
                        map(
                            lambda offset: str(offset + 1),
                            sorted(dnaMatch[key]["noCoverageOffsets"]),
                        )
                    ),
                )
            )
        ambiguousCount = len(dnaMatch[key]["ambiguousOffsets"])
        append(countPrint("%s    Ambiguous" % indent, ambiguousCount, length))
        extraCount = dnaMatch[key]["extraCount"]
        if extraCount:
            append(
                countPrint(
                    "%s    Extra nucleotides at end" % indent, extraCount, length
                )
            )

    if includeAmbiguousMatches and match["ambiguousMatches"]:
        append(f"{indent}Ambiguous matches:")
        for offset, a, b in match["ambiguousMatches"]:
            append(f"{indent}    {offset + 1} {a} {b}")

    if includeNonGapMismatches and match["nonGapMismatches"]:
        append(f"{indent}Non-gap mismatches:")
        for offset, a, b in match["nonGapMismatches"]:
            append(f"{indent}    {offset + 1} {a} {b}")

    return "\n".join(result)


def compareDNAReads(
    read1, read2, matchAmbiguous=True, gapChars="-", noCoverageChars=None, offsets=None
):
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
    @param noCoverageChars: An object supporting __contains__ with characters
        that indicate sequence positions that did not have coverage (for
        comparing sequences that may be consensuses from sequencing reads
        mapped to a reference.
    @param offsets: If not C{None}, a C{set} of offsets of interest. Offsets
        not in the set will not be considered.
    @return: A C{dict} with information about the match and the individual
        sequences (see below).
    """
    identicalMatchCount = ambiguousMatchCount = 0
    gapMismatchCount = nonGapMismatchCount = gapGapMismatchCount = 0
    noCoverageCount = noCoverageNoCoverageCount = 0
    read1ExtraCount = read2ExtraCount = 0
    read1GapOffsets = []
    read2GapOffsets = []
    read1AmbiguousOffsets = []
    read2AmbiguousOffsets = []
    read1NoCoverageOffsets = []
    read2NoCoverageOffsets = []
    empty = set()
    noCoverageChars = noCoverageChars or empty
    nonGapMismatches = []
    ambiguousMatches = []

    def _identicalMatch(a, b):
        """
        Two nucleotides are considered identical if they are the same character
        and they are not ambiguous.
        """
        return a == b and AMBIGUOUS.get(a, empty) == {a}

    def _ambiguousMatch(a, b, matchAmbiguous):
        """
        Checks if two characters match ambiguously if matchAmbiguous is True.
        A match is an ambiguous match if it is not an identical match, but the
        sets of ambiguous characters overlap.
        """
        return (
            matchAmbiguous
            and not _identicalMatch(a, b)
            and AMBIGUOUS.get(a, empty) & AMBIGUOUS.get(b, empty)
        )

    for offset, (a, b) in enumerate(
        zip_longest(read1.sequence.upper(), read2.sequence.upper())
    ):
        # Use 'is not None' in the following to allow an empty offsets set
        # to be passed.
        if offsets is not None and offset not in offsets:
            continue
        if len(AMBIGUOUS.get(a, "")) > 1:
            read1AmbiguousOffsets.append(offset)
        if len(AMBIGUOUS.get(b, "")) > 1:
            read2AmbiguousOffsets.append(offset)
        aNoCoverage = a in noCoverageChars
        bNoCoverage = b in noCoverageChars
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
        elif aNoCoverage or bNoCoverage:
            if aNoCoverage and bNoCoverage:
                noCoverageNoCoverageCount += 1
                read1NoCoverageOffsets.append(offset)
                read2NoCoverageOffsets.append(offset)
            elif aNoCoverage:
                noCoverageCount += 1
                read1NoCoverageOffsets.append(offset)
            else:
                noCoverageCount += 1
                read2NoCoverageOffsets.append(offset)
        else:
            # We have a character from both sequences (they could still be
            # gap characters).
            if a in gapChars:
                read1GapOffsets.append(offset)
                if b in gapChars:
                    # Both are gaps. This could happen if the sequences are
                    # taken from a multiple sequence alignment and some
                    # other sequence has resulted in a gap being inserted
                    # in both our seqeunces. This should never happen if
                    # our sequences were pairwise aligned (no sensible
                    # alignment program would have a reason to put a gap at
                    # the same place when aligning just two sequences).
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
                        ambiguousMatches.append((offset, a, b))
                    else:
                        nonGapMismatchCount += 1
                        nonGapMismatches.append((offset, a, b))

    return {
        "match": {
            "identicalMatchCount": identicalMatchCount,
            "ambiguousMatchCount": ambiguousMatchCount,
            "gapMismatchCount": gapMismatchCount,
            "gapGapMismatchCount": gapGapMismatchCount,
            "nonGapMismatchCount": nonGapMismatchCount,
            "noCoverageCount": noCoverageCount,
            "noCoverageNoCoverageCount": noCoverageNoCoverageCount,
            "ambiguousMatches": ambiguousMatches,
            "nonGapMismatches": nonGapMismatches,
        },
        "read1": {
            "ambiguousOffsets": read1AmbiguousOffsets,
            "extraCount": read1ExtraCount,
            "gapOffsets": read1GapOffsets,
            "noCoverageOffsets": read1NoCoverageOffsets,
        },
        "read2": {
            "ambiguousOffsets": read2AmbiguousOffsets,
            "extraCount": read2ExtraCount,
            "gapOffsets": read2GapOffsets,
            "noCoverageOffsets": read2NoCoverageOffsets,
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
            triplet = readSeq[offset : offset + 3]
            if triplet == "ATG":
                if readSeq[offset + 3] == "G":
                    if readSeq[offset - 3] in "GA":
                        kozakQualityCount = sum(
                            (
                                readSeq[offset - 1] == "C",
                                readSeq[offset - 2] == "C",
                                readSeq[offset - 4] == "C",
                                readSeq[offset - 5] == "C",
                                readSeq[offset - 6] == "G",
                            )
                        )

                        kozakQualityPercent = kozakQualityCount / 5.0 * 100
                        yield DNAKozakRead(
                            read, offset - 6, offset + 4, kozakQualityPercent
                        )
            offset += 1


class FloatBaseCounts:
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

        default = self._default = set("ACGT") if unknownAreAmbiguous else {"-"}

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
        fmt = "%d" if all(c == int(c) for b, c in self._sorted) else "%.2f"
        return "%s (%.3f)" % (
            " ".join(("%s:" + fmt) % (b, c) for b, c in self._sorted),
            self.highestFrequency(),
        )

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


def sequenceToRegex(sequence, wildcards="-?"):
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
            possible = "".join(sorted(AMBIGUOUS[base]))
        except KeyError:
            if base in wildcards:
                possible = "ACGT"
            else:
                raise

        append(("[%s]" % possible) if len(possible) > 1 else possible)

    return "".join(result)


def leastAmbiguous(bases):
    """
    Get the least ambiguous code for a set of DNA bases.

    @param bases: An iterable of C{str} DNA bases.
    @raise KeyError: If C[bases} contains an unknown nucleotide.
    @return: a C{str} DNA code (one of the values of the
        BASES_TO_AMBIGUOUS dict defined at top).
    """
    return BASES_TO_AMBIGUOUS["".join(sorted(set(map(str.upper, bases))))]


def leastAmbiguousFromCounts(bases, threshold):
    """
    Get the least ambiguous code given frequency counts for a set of DNA bases
    and a homogeneity threshold.

    @param bases: A C{dict} (or C{Counter}) mapping C{str} nucleotide
        letters to C{int} frequency counts.
    @param threshold: A C{float} homogeneity frequency.
    @raise ValueError: If any count or the thrshold is less than zero.
    @return: a C{str} DNA code (one of the values of the
        BASES_TO_AMBIGUOUS dict defined at top).
    """
    total = 0
    for base, count in bases.items():
        if count < 0:
            raise ValueError(f"Count for base {base!r} is negative ({count}).")
        total += count

    if threshold < 0:
        raise ValueError(f"Threshold cannot be negative ({threshold}).")

    if total == 0:
        return leastAmbiguous("ACGT")

    counts = sorted(bases.items(), key=itemgetter(1), reverse=True)
    cumulative = 0
    resultBases = set()
    for index in range(len(counts)):
        base, count = counts[index]
        cumulative += count
        resultBases.add(base)
        if cumulative / total >= threshold:
            break

    # Add bases with counts that tie the most-recently added base.
    lastCount = count
    for index in range(index + 1, len(counts)):
        base, count = counts[index]
        if count == lastCount:
            resultBases.add(base)
        else:
            break

    return leastAmbiguous(resultBases)


class Bases:
    """
    Manage a collection of (base, quality) pairs for a genome site.
    """

    __slots__ = ("count", "counts")

    def __init__(self):
        self.count = 0
        self.counts = dict.fromkeys("ACGT", 0)

    def __eq__(self, other):
        return self.count == other.count and self.counts == other.counts

    def __str__(self):
        return f"<Bases count={self.count}, bases={self.counts}"

    __repr__ = __str__

    def __getitem__(self, key):
        return self.counts[key]

    def __len__(self):
        return self.count

    def __add__(self, other):
        new = Bases()
        new.count = self.count + other.count
        for counts in self.counts, other.counts:
            for key, count in counts.items():
                new.counts[key] += count
        return new

    def __iadd__(self, other):
        self.count += other.count
        for key, count in other.counts.items():
            self.counts[key] += count
        return self

    def append(self, base, quality):
        """
        Append a (base, quality) pair.

        @param base: A C{str} nucleotide base.
        @param quality: An C{int} nucleotide quality score.
        @return: C{self}.
        """
        if base != "N":
            self.count += 1
            self.counts[base] += quality

        return self

    def consensus(self, threshold, minCoverage, lowCoverage, noCoverage):
        """
        Get the base that can be used as part of a consensus.

        If there are sufficient reads, this is the least-ambiguous nucleotide
        code for our bases, given a required homogeneity threshold. Otherwise,
        the low coverage value.

        @param threshold: A C{float} threshold. This fraction, at least, of the
            most-common nucleotides at a site are used to determine the
            consensus nucleotide (or ambiguous symbol if more than one
            nucleotide is required to achieve this threshold). If there is a
            tie in nucleotide counts at a site that causes the threshold to be
            met, all nucleotides of equeal frequncy will be included in the
            ambiguous symbol for that site. This is perhaps better explained
            with an example. See
            https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
            and the corresponding testGeneiousExamplesTie test in
            test/test_dna.py
        @param minCoverage: An C{int} minimum number of reads that must cover a
            site for a consensus base to be called. If fewer reads cover a
            site, the C{lowCoverage} value is used.
        @param lowCoverage: A C{str} indicating what base to use when
            0 < N < minCoverage reads cover a site.
        @param noCoverage: A C{str} indicating what base to use when
            no reads cover the site.
        @return: A C{str} nucleotide code. This will be an ambiguous code if
            the homogeneity C{threshold} is not met.
        """
        return (
            noCoverage
            if self.count == 0
            else (
                lowCoverage
                if self.count < minCoverage
                else leastAmbiguousFromCounts(self.counts, threshold)
            )
        )

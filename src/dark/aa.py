from typing import Generator, Union, Iterable, Mapping, Optional
from itertools import zip_longest

from dark.aaVars import (
    ABBREV3,
    ABBREV3_TO_ABBREV1,
    CODONS,
    NAMES,
    NAMES_TO_ABBREV1,
    PROPERTIES,
    PROPERTY_CLUSTERS,
    PROPERTY_DETAILS,
)

from dark.reads import Read
from dark.utils import countPrint


class AminoAcid:
    """
    Hold information about an amino acid.

    @param name: The full C{str} name of the amino acid.
    @param abbrev3: The 3-letter C{str} abbreviation of the amino acid,
        e.g., 'Arg'.
    @param abbrev1: The 1-letter C{str} abbreviation of the amino acid,
        e.g., 'A'.
    @param codons: A C{list} of 3-letter codons for the amino acid.
    @param properties: An C{int} logical-AND of the various properties
        (see PROPERTIES, above) of this amino acid.
    @param propertyDetails: A C{dict} containing property names and values
        for this amino acid. E.g.:
        {
            'aliphaticity': -0.157024793388,
            'aromaticity': -0.0642673521851,
            'composition': -0.527272727273,
            'hydrogenation': -0.401797175866,
            'hydropathy': -1.0,
            'hydroxythiolation': -0.51486325802,
            'iep': 1.0,
            'polar requirement': 0.0487804878049,
            'polarity': 0.382716049383,
            'volume': 0.449101796407,
        }
    @param propertyClusters: A C{dict} containing the property names and
        clusters for this amino acid. E.g.:
        {
            'aliphaticity': 1,
            'aromaticity': 2,
            'composition': 1,
            'hydrogenation': 1,
            'hydropathy': 2,
            'hydroxythiolation': 4,
            'iep': 2,
            'polar requirement': 1,
            'polarity': 1,
            'volume': 4,
        }
    """

    def __init__(
        self,
        name: str,
        abbrev3: str,
        abbrev1: str,
        codons: tuple[str, ...],
        properties: int,
        propertyDetails: dict[str, float],
        propertyClusters: dict[str, int],
    ):
        self.name = name
        self.abbrev3 = abbrev3
        self.abbrev1 = abbrev1
        self.codons = codons
        self.properties = properties
        self.propertyDetails = propertyDetails
        self.propertyClusters = propertyClusters


def find(s: str) -> Generator[AminoAcid, None, None]:
    """
    Find an amino acid whose name or abbreviation is s.

    @param s: A C{str} amino acid specifier. This may be a full name,
        a 3-letter abbreviation or a 1-letter abbreviation. Case is ignored.
    return: A generator that yields matching L{AminoAcid} instances.
    """

    abbrev1 = None
    origS = s

    if " " in s:
        # Convert first word to title case, others to lower.
        first, rest = s.split(" ", 1)
        s = first.title() + " " + rest.lower()
    else:
        s = s.title()

    if s in NAMES:
        abbrev1 = s
    elif s in ABBREV3_TO_ABBREV1:
        abbrev1 = ABBREV3_TO_ABBREV1[s]
    elif s in NAMES_TO_ABBREV1:
        abbrev1 = NAMES_TO_ABBREV1[s]
    else:
        # Look for a 3-letter codon.
        def findCodon(target):
            for abbrev1, codons in CODONS.items():
                for codon in codons:
                    if codon == target:
                        return abbrev1

        # Allow for RNA in this lookup by replacing U with T.
        abbrev1 = findCodon(origS.upper().replace("U", "T"))

    if abbrev1:
        abbrev1s = [abbrev1]
    else:
        # Try partial matching on names.
        abbrev1s = []
        sLower = s.lower()
        for abbrev1, name in NAMES.items():
            if name.lower().find(sLower) > -1:
                abbrev1s.append(abbrev1)

    for abbrev1 in abbrev1s:
        yield AminoAcid(
            NAMES[abbrev1],
            ABBREV3[abbrev1],
            abbrev1,
            CODONS[abbrev1],
            PROPERTIES[abbrev1],
            PROPERTY_DETAILS[abbrev1],
            PROPERTY_CLUSTERS[abbrev1],
        )


def _propertiesOrClustersForSequence(
    sequence: Read,
    propertyNames: Iterable[str],
    propertyValues: Mapping[str, Mapping[str, Union[int, float]]],
    missingAAValue: Union[int, float],
) -> dict[str, list[Union[int, float]]]:
    """
    Extract amino acid property values or cluster numbers for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{propertyValues} C{dict}).
    @param propertyValues: A C{dict} in the form of C{PROPERTY_DETAILS} or
        C{PROPERTY_CLUSTERS} (see above).
    @param missingAAValue: A C{float} value to use for properties when an AA
        (e.g., 'X') is not known.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: A C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding C{float} property value in C{propertyValues} in
        order of sequence position.
    """
    propertyNames = sorted(map(str.lower, set(propertyNames)))

    # Make sure all mentioned property names exist for at least one AA.
    knownProperties: set[str] = set()
    for names in propertyValues.values():
        knownProperties.update(names)
    unknown = set(propertyNames) - knownProperties
    if unknown:
        raise ValueError(
            "Unknown propert%s: %s."
            % ("y" if len(unknown) == 1 else "ies", ", ".join(unknown))
        )

    aas = sequence.sequence.upper()
    result: dict[str, list[float]] = {}

    for propertyName in propertyNames:
        result[propertyName] = []
        append = result[propertyName].append
        for aa in aas:
            try:
                properties = propertyValues[aa]
            except KeyError:
                # No such AA.
                append(missingAAValue)
            else:
                append(properties[propertyName])

    return result


def propertiesForSequence(
    sequence: Read, propertyNames: Iterable[str], missingAAValue: float = -1.1
):
    """
    Extract amino acid property values for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{dark.aa.PROPERTY_DETAILS} C{dict}).
    @param missingAAValue: A C{float} value to use for properties when an AA
        (e.g., 'X') is not known.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: A C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding property value according to sequence
        position.
    """
    return _propertiesOrClustersForSequence(
        sequence, propertyNames, PROPERTY_DETAILS, missingAAValue
    )


def clustersForSequence(
    sequence: Read, propertyNames: Iterable[str], missingAAValue: int = 0
):
    """
    Extract amino acid property cluster numbers for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{dark.aa.PROPERTY_CLUSTERS} C{dict}).
    @param missingAAValue: An C{int} value to use for properties when an AA
        (e.g., 'X') is not known.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: A C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding C{int} property cluster number according to
        sequence position.
    """
    return _propertiesOrClustersForSequence(
        sequence, propertyNames, PROPERTY_CLUSTERS, missingAAValue
    )


# It would be nice, in some universe, to do the following. But we access the
# match dict using a variable key and mypy cannot check that, so this doesn't work.
# Leaving it here in case it can be used for stricter checking in the future.
#
# _OverallMatch = TypedDict(
#     "_OverallMatch",
#     {
#         "matchCount": int,
#         "gapMismatchCount": int,
#         "gapGapMismatchCount": int,
#         "nonGapMismatchCount": int,
#     },
# )
#
# _Read1Match = TypedDict(
#     "_Read1Match",
#     {
#         "extraCount": int,
#         "gapOffsets": list[int],
#     },
# )
#
# _Read2Match = TypedDict(
#     "_Read2Match",
#     {
#         "extraCount": int,
#         "gapOffsets": list[int],
#     },
# )
#
# Match = TypedDict(
#     "Match",
#     {
#         "match": _OverallMatch,
#         "read1": _Read1Match,
#         "read2": _Read2Match,
#     },
# )


# Be more relaxed.
_Match = Mapping[str, Mapping]


def matchToString(
    aaMatch: _Match,
    read1: Read,
    read2: Read,
    indent: str = "",
    offsets: Optional[set[int]] = None,
) -> str:
    """
    Format amino acid sequence match as a string.

    @param aaMatch: A C{dict} returned by C{compareAaReads}.
    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param indent: A C{str} to indent all returned lines with.
    @param offsets: If not C{None}, a C{set} of offsets of interest that were
        only considered when making C{match}.
    @return: A C{str} describing the match.
    """
    match = aaMatch["match"]
    matchCount = match["matchCount"]
    gapMismatchCount = match["gapMismatchCount"]
    gapGapMismatchCount = match["gapGapMismatchCount"]
    nonGapMismatchCount = match["nonGapMismatchCount"]

    if offsets:
        len1 = len2 = len(offsets)
    else:
        len1, len2 = map(len, (read1, read2))

    result: list[str] = []
    append = result.append

    append("%sMatches:" % indent)
    append(countPrint("%s  Overall" % indent, matchCount, len1, len2))
    append(
        countPrint(
            "%s  Not involving gaps (i.e., identities)" % indent,
            matchCount,
            matchCount + nonGapMismatchCount,
        )
    )

    append("%sMismatches:" % indent)
    mismatchCount = gapMismatchCount + gapGapMismatchCount + nonGapMismatchCount
    append(countPrint("%s  Overall" % indent, mismatchCount, len1, len2))
    append(
        countPrint(
            "%s  Not involving gaps (i.e., conflicts)" % indent,
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

    for read, key in zip((read1, read2), ("read1", "read2")):
        append("%sId: %s" % (indent, read.id))
        length = len(read)
        append("%s  Length: %d" % (indent, length))
        gapCount = len(aaMatch[key]["gapOffsets"])
        append(countPrint("%s  Gaps" % indent, gapCount, length))
        if gapCount:
            append(
                "%s  Gap locations (1-based): %s"
                % (
                    indent,
                    ", ".join(
                        map(
                            lambda offset: str(offset + 1),
                            sorted(aaMatch[key]["gapOffsets"]),
                        )
                    ),
                )
            )
        extraCount = aaMatch[key]["extraCount"]
        if extraCount:
            append(
                countPrint("%s  Extra amino acids at end" % indent, extraCount, length)
            )

    return "\n".join(result)


def compareAaReads(
    read1: Read, read2: Read, gapChars: str = "-", offsets: Optional[set[int]] = None
) -> _Match:
    """
    Compare two amino acid sequences.

    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param gapChars: An object supporting __contains__ with characters that
        should be considered to be gaps.
    @param offsets: If not C{None}, a C{set} of offsets of interest. Offsets
        not in the set will not be considered.
    @return: A C{dict} with information about the match and the individual
        sequences (see below).
    """
    matchCount = 0
    gapMismatchCount = nonGapMismatchCount = gapGapMismatchCount = 0
    read1ExtraCount = read2ExtraCount = 0
    read1GapOffsets = []
    read2GapOffsets = []

    for offset, (a, b) in enumerate(
        zip_longest(read1.sequence.upper(), read2.sequence.upper())
    ):
        # Use 'is not None' in the following to allow an empty offsets set
        # to be passed.
        if offsets is not None and offset not in offsets:
            continue
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
                    if a == b:
                        matchCount += 1
                    else:
                        nonGapMismatchCount += 1

    return {
        "match": {
            "matchCount": matchCount,
            "gapMismatchCount": gapMismatchCount,
            "gapGapMismatchCount": gapGapMismatchCount,
            "nonGapMismatchCount": nonGapMismatchCount,
        },
        "read1": {
            "extraCount": read1ExtraCount,
            "gapOffsets": read1GapOffsets,
        },
        "read2": {
            "extraCount": read2ExtraCount,
            "gapOffsets": read2GapOffsets,
        },
    }

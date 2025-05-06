from __future__ import annotations

import re
from math import ceil
from typing import Optional, TYPE_CHECKING
import argparse

if TYPE_CHECKING:
    from dark.reads import Reads
    from dark.titles import TitleAlignments
from dark.simplify import simplifyTitle
from dark.utils import parseRangeExpression


class TitleFilter:
    """
    Provide an acceptance test for sequence titles.

    @param whitelist: If not C{None}, a C{set} of exact titles that are always
        acceptable.
    @param blacklist: If not C{None}, a C{set} of exact titles that are never
        acceptable.
    @param whitelistFile: If not C{None}, a C{str} filename containing lines
        that give exact titles that are always acceptable.
    @param blacklistFile: If not C{None}, a C{str} filename containing lines
        that give exact titles that are never acceptable.
    @param positiveRegex: If not C{None}, a C{str} regex that sequence titles
        must match (case is ignored).
    @param negativeRegex: If not C{None}, a C{str} regex that sequence titles
        must not match (case is ignored).
    @param truncateAfter: A C{str} that titles will be truncated beyond. If
        a truncated title has already been seen, that title will no longer
        be acceptable.
    """

    # The following are return values for the self.accept method.
    REJECT = 0
    WHITELIST_ACCEPT = 1
    DEFAULT_ACCEPT = 2

    def __init__(
        self,
        whitelist: Optional[set[str]] = None,
        blacklist: Optional[set[str]] = None,
        whitelistFile: Optional[str] = None,
        blacklistFile: Optional[str] = None,
        positiveRegex: Optional[str] = None,
        negativeRegex: Optional[str] = None,
        truncateAfter: Optional[str] = None,
    ) -> None:
        whitelist = whitelist or set()
        if whitelistFile:
            with open(whitelistFile) as fp:
                for line in fp:
                    whitelist.add(line.strip())
        self._whitelist = whitelist

        blacklist = blacklist or set()
        if blacklistFile:
            with open(blacklistFile) as fp:
                for line in fp:
                    blacklist.add(line.strip())
        self._blacklist = blacklist

        if truncateAfter is None:
            self._truncated: Optional[dict[str, str]] = None
        else:
            self._truncateAfter = truncateAfter
            self._truncated = {}

        if positiveRegex is None:
            self._positiveRegex = None
        else:
            self._positiveRegex = re.compile(positiveRegex, re.I)

        if negativeRegex is None:
            self._negativeRegex = None
        else:
            self._negativeRegex = re.compile(negativeRegex, re.I)

    def accept(self, title: str) -> int:
        """
        Return a value (see below) to indicate if a title is acceptable (and,
        if so, in what way).

        @param title: A C{str} sequence title.
        @return: An C{int} to indicate an acceptable title or not. This will be

            C{self.REJECT} if the title is unacceptable.
            C{self.WHITELIST_ACCEPT} if the title is whitelisted.
            C{self.DEFAULT_ACCEPT} if the title is acceptable by default.

            These three values are needed so our caller can distinguish between
            the two reasons for acceptance.
        """
        if self._whitelist and title in self._whitelist:
            return self.WHITELIST_ACCEPT

        if self._blacklist and title in self._blacklist:
            return self.REJECT

        # If we have a positive regex but we don't match it, reject.
        if self._positiveRegex and self._positiveRegex.search(title) is None:
            return self.REJECT

        # If we have a negative regex and we do match it, reject.
        if self._negativeRegex and self._negativeRegex.search(title) is not None:
            return self.REJECT

        if self._truncated is not None:
            truncated = simplifyTitle(title, self._truncateAfter)
            if truncated in self._truncated:
                # We've already seen this (truncated) title. Reject unless
                # this is the original title that we truncated to make this
                # entry. That title must continue to be accepted.
                if self._truncated[truncated] == title:
                    return self.DEFAULT_ACCEPT
                else:
                    return self.REJECT
            else:
                self._truncated[truncated] = title

        return self.DEFAULT_ACCEPT


class SequenceFilter:
    """
    Provide an acceptance test for sequences.

    @param whitelist: If not C{None}, a C{set} of exact sequences that are always
        acceptable.
    @param blacklist: If not C{None}, a C{set} of exact sequences that are never
        acceptable.
    @param whitelistFile: If not C{None}, a C{str} filename containing lines
        that give exact sequences that are always acceptable.
    @param blacklistFile: If not C{None}, a C{str} filename containing lines
        that give exact sequences that are never acceptable.
    @param positiveRegex: If not C{None}, a C{str} regex that sequences
        must match (case is ignored).
    @param negativeRegex: If not C{None}, a C{str} regex that sequences
        must not match (case is ignored).
    """

    # The following are return values for the self.accept method.
    REJECT = 0
    WHITELIST_ACCEPT = 1
    DEFAULT_ACCEPT = 2

    def __init__(
        self,
        whitelist: Optional[set[str]] = None,
        blacklist: Optional[set[str]] = None,
        whitelistFile: Optional[str] = None,
        blacklistFile: Optional[str] = None,
        positiveRegex: Optional[str] = None,
        negativeRegex: Optional[str] = None,
    ) -> None:
        whitelist = whitelist or set()
        if whitelistFile:
            with open(whitelistFile) as fp:
                for line in fp:
                    whitelist.add(line.strip())
        self._whitelist = whitelist

        blacklist = blacklist or set()
        if blacklistFile:
            with open(blacklistFile) as fp:
                for line in fp:
                    blacklist.add(line.strip())
        self._blacklist = blacklist

        if positiveRegex is None:
            self._positiveRegex = None
        else:
            self._positiveRegex = re.compile(positiveRegex, re.I)

        if negativeRegex is None:
            self._negativeRegex = None
        else:
            self._negativeRegex = re.compile(negativeRegex, re.I)

    def accept(self, sequence: str) -> int:
        """
        Return a value (see below) to indicate if a sequence is acceptable (and,
        if so, in what way).

        @param sequence: A C{str} sequence.
        @return: An C{int} to indicate an acceptable sequence or not. This will be

            C{self.REJECT} if the sequence is unacceptable.
            C{self.WHITELIST_ACCEPT} if the sequence is whitelisted.
            C{self.DEFAULT_ACCEPT} if the sequence is acceptable by default.

            These three values are needed so our caller can distinguish between
            the two reasons for acceptance.
        """
        if self._whitelist and sequence in self._whitelist:
            return self.WHITELIST_ACCEPT

        if self._blacklist and sequence in self._blacklist:
            return self.REJECT

        # If we have a positive regex but we don't match it, reject.
        if self._positiveRegex and self._positiveRegex.search(sequence) is None:
            return self.REJECT

        # If we have a negative regex and we do match it, reject.
        if self._negativeRegex and self._negativeRegex.search(sequence) is not None:
            return self.REJECT

        return self.DEFAULT_ACCEPT


class ReadSetFilter:
    """
    Provide an acceptance test based on sequence read set.

    @param minNew: The C{float} fraction of its reads by which a new read set
        must differ from all previously seen read sets in order to be
        considered acceptably different.
    """

    def __init__(self, minNew: float) -> None:
        self._minNew = minNew
        self._titles: dict[str, tuple[set[str], list[str]]] = {}

    def accept(self, title: str, titleAlignments: TitleAlignments) -> bool:
        """
        Return C{True} if the read id set in C{titleAlignments} is sufficiently
        different from all previously seen read sets.

        @param title: A C{str} sequence title.
        @param titleAlignments: An instance of L{TitleAlignment}.
        @return: A C{bool} indicating whether a title has an acceptably novel
            read set or not.
        """

        # Sanity check: titles can only be passed once.
        assert title not in self._titles, "Title %r seen multiple times." % title

        readIds = titleAlignments.readIds()
        newReadsRequired = ceil(self._minNew * len(readIds))

        for readSet, invalidatedTitles in self._titles.values():
            if len(readIds - readSet) < newReadsRequired:
                # Add this title to the set of titles invalidated by this
                # previously seen read set.
                invalidatedTitles.append(title)
                return False

        # Remember the new read set and an empty list of invalidated titles.
        self._titles[title] = (readIds, [])

        return True

    def invalidates(self, title: str) -> list[str]:
        """
        Report on which other titles were invalidated by a given title.

        @param title: A C{str} sequence title.
        @return: A C{list} of titles that the passed title invalidated.
        """
        try:
            return self._titles[title][1]
        except KeyError:
            return []


def addFASTAFilteringCommandLineOptions(parser: argparse.ArgumentParser) -> None:
    """
    Add standard FASTA filtering command-line options to an argparse parser.

    These are options that can be used to select or omit entire FASTA records,
    NOT options that change them (for that see
    addFASTAEditingCommandLineOptions).

    @param parser: An C{argparse.ArgumentParser} instance.
    """
    parser.add_argument(
        "--minLength", type=int, metavar="N", help="The minimum sequence length."
    )

    parser.add_argument(
        "--maxLength", type=int, metavar="N", help="The maximum sequence length."
    )

    parser.add_argument(
        "--maxNFraction",
        type=float,
        metavar="N",
        help="The maximum fraction of Ns that can be present in the sequence.",
    )

    parser.add_argument(
        "--whitelist",
        action="append",
        metavar="SEQUENCE-ID",
        help="Sequence titles (ids) that should be whitelisted. May be repeated.",
    )

    parser.add_argument(
        "--blacklist",
        action="append",
        metavar="SEQUENCE-ID",
        help="Sequence titles (ids) that should be blacklisted. May be repeated.",
    )

    parser.add_argument(
        "--whitelistFile",
        metavar="SEQUENCE-ID-FILE",
        help=(
            "The name of a file that contains sequence titles (ids) that "
            "should be whitelisted, one per line"
        ),
    )

    parser.add_argument(
        "--blacklistFile",
        metavar="SEQUENCE-ID-FILE",
        help=(
            "The name of a file that contains sequence titles (ids) that "
            "should be blacklisted, one per line."
        ),
    )

    parser.add_argument(
        "--titleRegex",
        metavar="REGEX",
        help="A regex that sequence titles (ids) must match.",
    )

    parser.add_argument(
        "--negativeTitleRegex",
        metavar="REGEX",
        help="A regex that sequence titles (ids) must not match.",
    )

    parser.add_argument(
        "--sequenceWhitelist",
        action="append",
        metavar="SEQUENCE",
        help="Sequences that should be whitelisted. May be repeated.",
    )

    parser.add_argument(
        "--sequenceBlacklist",
        action="append",
        metavar="SEQUENCE",
        help="Sequences that should be blacklisted. May be repeated.",
    )

    parser.add_argument(
        "--sequenceWhitelistFile",
        metavar="SEQUENCE-FILE",
        help=(
            "The name of a file that contains sequences that should be whitelisted, "
            "one per line."
        ),
    )

    parser.add_argument(
        "--sequenceBlacklistFile",
        metavar="SEQUENCE-FILE",
        help=(
            "The name of a file that contains sequences that should be blacklisted, "
            "one per line."
        ),
    )

    parser.add_argument(
        "--sequenceRegex",
        metavar="REGEX",
        help="A regex that sequences must match.",
    )

    parser.add_argument(
        "--sequenceNegativeRegex",
        metavar="REGEX",
        help="A regex that sequences must not match.",
    )

    # A mutually exclusive group for --keepSequences and --removeSequences.
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--keepSequences",
        metavar="NUMBER,RANGE,...",
        help=(
            "Specify (1-based) ranges of sequence numbers that should be "
            "kept. E.g., --keepSequences 1-3,5 will output just the 1st, "
            "2nd, 3rd, and 5th sequences. All others will be omitted."
        ),
    )

    group.add_argument(
        "--removeSequences",
        metavar="NUMBER,RANGE,...",
        help=(
            "Specify (1-based) ranges of sequence numbers that should be "
            "removed. E.g., --removeSequences 1-3,5 will output all but the "
            "1st, 2nd, 3rd, and 5th sequences. All others will be ouput."
        ),
    )

    parser.add_argument(
        "--head",
        type=int,
        metavar="N",
        help="Only the first N sequences will be printed.",
    )

    # A mutually exclusive group for --removeDuplicates and
    # --removeDuplicatesById
    _removeGroup = parser.add_mutually_exclusive_group()

    _removeGroup.add_argument(
        "--removeDuplicates",
        action="store_true",
        help=(
            "Duplicate reads will be removed, based only on sequence identity. The "
            "first occurrence is kept."
        ),
    )

    _removeGroup.add_argument(
        "--removeDuplicatesById",
        action="store_true",
        help=(
            "Duplicate reads will be removed, based only on read id. The first "
            "occurrence is kept."
        ),
    )

    parser.add_argument(
        "--removeDuplicatesUseMD5",
        action="store_true",
        help=(
            "MD5 sums will be stored instead of the full sequence or read "
            "id when either --removeDuplicates or removeDuplicatesById are "
            "given. One of --removeDuplicates and --removeDuplicatesById "
            "must also be given. Note that this makes duplicate removal "
            "probabilistic. This option can be used to reduce the amount "
            "of RAM consumed during duplicate removal."
        ),
    )

    # See the docstring for dark.reads.Reads.filter for more detail on
    # randomSubset.
    parser.add_argument(
        "--randomSubset",
        type=int,
        metavar="N",
        help=(
            "An integer giving the number of sequences that should be kept. "
            "These will be selected at random."
        ),
    )

    # See the docstring for dark.reads.Reads.filter for more detail on
    # trueLength.
    parser.add_argument(
        "--trueLength",
        type=int,
        metavar="N",
        help=(
            "The number of reads in the FASTA input. Only to be used with "
            "randomSubset"
        ),
    )

    parser.add_argument(
        "--sampleFraction",
        type=float,
        metavar="FRACTION",
        help=(
            "A [0.0, 1.0] C{float} indicating a fraction of the reads that "
            "should be allowed to pass through the filter. The sample size "
            "will only be approximately the product of the sample fraction "
            "and the number of reads. The sample is taken at random."
        ),
    )

    parser.add_argument(
        "--sequenceNumbersFile",
        metavar="FILENAME",
        help=(
            "A file of (1-based) sequence numbers to retain. Numbers must "
            "be one per line."
        ),
    )


def parseFASTAFilteringCommandLineOptions(
    args: argparse.Namespace, reads: Reads
) -> Reads:
    """
    Examine parsed FASTA filtering command-line options and return filtered
    reads.

    @param args: An argparse namespace, as returned by the argparse
        C{parse_args} function.
    @param reads: A C{Reads} instance to filter.
    @return: The filtered C{Reads} instance.
    """
    keepSequences = (
        parseRangeExpression(args.keepSequences, convertToZeroBased=True)
        if args.keepSequences
        else None
    )

    removeSequences = (
        parseRangeExpression(args.removeSequences, convertToZeroBased=True)
        if args.removeSequences
        else None
    )

    return reads.filter(
        minLength=args.minLength,
        maxLength=args.maxLength,
        maxNFraction=args.maxNFraction,
        whitelist=set(args.whitelist) if args.whitelist else None,
        blacklist=set(args.blacklist) if args.blacklist else None,
        whitelistFile=args.whitelistFile,
        blacklistFile=args.blacklistFile,
        titleRegex=args.titleRegex,
        negativeTitleRegex=args.negativeTitleRegex,
        sequenceWhitelist=set(args.sequenceWhitelist)
        if args.sequenceWhitelist
        else None,
        sequenceBlacklist=set(args.sequenceBlacklist)
        if args.sequenceBlacklist
        else None,
        sequenceWhitelistFile=args.sequenceWhitelistFile,
        sequenceBlacklistFile=args.sequenceBlacklistFile,
        sequenceRegex=args.sequenceRegex,
        sequenceNegativeRegex=args.sequenceNegativeRegex,
        keepSequences=keepSequences,
        removeSequences=removeSequences,
        head=args.head,
        removeDuplicates=args.removeDuplicates,
        removeDuplicatesById=args.removeDuplicatesById,
        removeDuplicatesUseMD5=args.removeDuplicatesUseMD5,
        randomSubset=args.randomSubset,
        trueLength=args.trueLength,
        sampleFraction=args.sampleFraction,
        sequenceNumbersFile=args.sequenceNumbersFile,
    )


def addFASTAEditingCommandLineOptions(parser: argparse.ArgumentParser) -> None:
    """
    Add standard FASTA editing command-line options to an argparse parser.

    These are options that can be used to alter FASTA records, NOT options
    that simply select or reject those things (for those see
    addFASTAFilteringCommandLineOptions).

    @param parser: An C{argparse.ArgumentParser} instance.
    """
    # A mutually exclusive group for --keepSites, --keepSitesFile,
    # --removeSites, and --removeSitesFile.
    group = parser.add_mutually_exclusive_group()

    # In the 4 options below, the 'indices' alternate names are kept for
    # backwards compatibility.
    group.add_argument(
        "--keepSites",
        "--keepIndices",
        help=(
            "Specify 1-based sequence sites to keep. All other sites will "
            "be removed. The sites must be given in the form e.g., "
            "24,100-200,260. Note that the requested sites will be taken "
            "from the input sequences in order, not in the order given by "
            "--keepSites. I.e., --keepSites 5,8-10 will get you the same "
            "result as --keepSites 8-10,5."
        ),
    )

    group.add_argument(
        "--keepSitesFile",
        "--keepIndicesFile",
        help=(
            "Specify a file containing 1-based sites to keep. All other "
            "sequence sites will be removed. Lines in the file must be "
            "given in the form e.g., 24,100-200,260. See --keepSites for "
            "more detail."
        ),
    )

    group.add_argument(
        "--removeSites",
        "--removeIndices",
        help=(
            "Specify 1-based sites to remove. All other sequence sites will "
            "be kept. The sites must be given in the form e.g., "
            "24,100-200,260. See --keepSites for more detail."
        ),
    )

    group.add_argument(
        "--removeSitesFile",
        "--removeIndicesFile",
        help=(
            "Specify a file containing 1-based sites to remove. All other "
            "sequence sites will be kept. Lines in the file must be given "
            "in the form e.g., 24,100-200,260. See --keepSites for more "
            "detail."
        ),
    )

    parser.add_argument(
        "--removeGaps",
        action="store_true",
        help="If True, gap ('-') characters in sequences will be removed.",
    )

    parser.add_argument(
        "--truncateTitlesAfter",
        help=(
            "A string that sequence titles (ids) will be truncated beyond. "
            "If the truncated version of a title has already been seen, "
            "that title will be skipped."
        ),
    )

    parser.add_argument(
        "--removeDescriptions",
        action="store_true",
        help=(
            "Read id descriptions will be removed. The "
            "description is the part of a sequence id after the "
            "first whitespace (if any)."
        ),
    )

    parser.add_argument(
        "--idLambda",
        metavar="LAMBDA-FUNCTION",
        help=(
            "A one-argument function taking and returning a read id. "
            "E.g., --idLambda \"lambda id: id.split('_')[0]\" or "
            '--idLambda "lambda id: id[:10]". If the function returns None, '
            "the read will be filtered out."
        ),
    )

    parser.add_argument(
        "--readLambda",
        metavar="LAMBDA-FUNCTION",
        help=(
            "A one-argument function taking and returning a read. "
            "E.g., --readLambda \"lambda r: Read(r.id.split('_')[0], "
            "r.sequence.strip('-'))\". Make sure to also modify the "
            "quality string if you change the length of a FASTQ sequence. "
            "If the function returns None, the read will be filtered out. "
            "The function will be passed to eval with the dark.reads "
            "classes Read, DNARead, AARead, etc. all in scope."
        ),
    )

    parser.add_argument(
        "--reverse",
        action="store_true",
        help="Reverse the sequences. Note that this is NOT reverse complementing.",
    )

    parser.add_argument(
        "--reverseComplement",
        action="store_true",
        help="Reverse complement the sequences.",
    )

    parser.add_argument(
        "--rotate",
        type=int,
        metavar="N",
        help=(
            "Rotate sequences. A positive integer rotates the sequence to the right, "
            "negative to the left"
        ),
    )

    seqCaseGroup = parser.add_mutually_exclusive_group()

    seqCaseGroup.add_argument(
        "--upper",
        action="store_true",
        help="Convert sequences to uppercase.",
    )

    seqCaseGroup.add_argument(
        "--lower",
        action="store_true",
        help="Convert sequences to lowercase.",
    )

    idCaseGroup = parser.add_mutually_exclusive_group()

    idCaseGroup.add_argument(
        "--upperId",
        action="store_true",
        help="Convert sequence IDs to uppercase.",
    )

    idCaseGroup.add_argument(
        "--lowerId",
        action="store_true",
        help="Convert sequence IDs to lowercase.",
    )


def parseFASTAEditingCommandLineOptions(
    args: argparse.Namespace, reads: Reads
) -> Reads:
    """
    Examine parsed FASTA editing command-line options and return information
    about kept sites and sequences.

    @param args: An argparse namespace, as returned by the argparse
        C{parse_args} function.
    @param reads: A C{Reads} instance to filter.
    @return: The filtered C{Reads} instance.
    """
    removeGaps = args.removeGaps
    removeDescriptions = args.removeDescriptions
    truncateTitlesAfter = args.truncateTitlesAfter
    keepSites = (
        parseRangeExpression(args.keepSites, convertToZeroBased=True)
        if args.keepSites
        else None
    )

    if args.keepSitesFile:
        keepSites = keepSites or set()
        with open(args.keepSitesFile) as fp:
            for lineNumber, line in enumerate(fp):
                try:
                    keepSites.update(
                        parseRangeExpression(line, convertToZeroBased=True)
                    )
                except ValueError as e:
                    raise ValueError(
                        "Keep sites file %r line %d could not be parsed: "
                        "%s" % (args.keepSitesFile, lineNumber, e)
                    )

    removeSites = (
        parseRangeExpression(args.removeSites, convertToZeroBased=True)
        if args.removeSites
        else None
    )

    if args.removeSitesFile:
        removeSites = removeSites or set()
        with open(args.removeSitesFile) as fp:
            for lineNumber, line in enumerate(fp):
                try:
                    removeSites.update(
                        parseRangeExpression(line, convertToZeroBased=True)
                    )
                except ValueError as e:
                    raise ValueError(
                        "Remove sites file %r line %d parse error: %s"
                        % (args.removeSitesFile, lineNumber, e)
                    )

    return reads.filter(
        removeGaps=removeGaps,
        truncateTitlesAfter=truncateTitlesAfter,
        removeDescriptions=removeDescriptions,
        idLambda=args.idLambda,
        readLambda=args.readLambda,
        keepSites=keepSites,
        removeSites=removeSites,
        reverse=args.reverse,
        reverseComplement=args.reverseComplement,
        upper=args.upper,
        lower=args.lower,
        upperId=args.upperId,
        lowerId=args.lowerId,
        rotate=args.rotate,
    )

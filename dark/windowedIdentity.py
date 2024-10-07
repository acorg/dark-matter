import re
from operator import attrgetter
from typing import Optional
from collections import defaultdict
from collections.abc import Iterable

from dark.dna import compareDNAReads
from dark.reads import DNARead


class WindowedIdentity:
    """
    Compute windowed nucleotide identity for a set of reads against a reference.
    """

    def __init__(self, reads: Iterable[DNARead]) -> None:
        self._reads = reads
        self._validated = False

    def _checkReads(self) -> None:
        """
        Check all reads have the same length.

        @param reads: A C{list} of C{DNARead} instances.
        """
        if self._validated:
            return

        if not self._reads:
            raise ValueError("Empty FASTA input!")

        first = True

        for count, read in enumerate(self._reads, start=1):
            if first:
                length = len(read)
                if length == 0:
                    raise ValueError("The first input sequence has length zero!")
                first = False
            else:
                if len(read) != length:
                    raise ValueError(
                        f"Sequence number {count}, with id {read.id!r}, has length "
                        f"{len(read)} which is not equal to the previous input sequence "
                        f"length(s) of {length}."
                    )

        self._validated = True

    def _matchingReads(self, pattern: str, reads: Iterable[DNARead]) -> list[DNARead]:
        """
        Find the reads whose ids match a pattern. The matching looks for the regular
        expression anywhere in the read ids.

        @param pattern: A regular expression.
        @param reads: The reads to examine.
        @return: A list of matching reads.
        """
        regex = re.compile(pattern)
        return [read for read in reads if regex.search(read.id)]

    def getIdentity(
        self,
        referenceRegex: str,
        window: int,
        minWindow: int,
        jump: int,
        includeRegex: str = ".",
        sort: bool = False,
        strict: bool = True,
        startOffset: int = 0,
        stopOffset: Optional[int] = None,
    ) -> tuple[DNARead, list[int], dict[DNARead, list[float]]]:
        """
        Compute windowed mucleotide identity of the reads whose ids are matched by
        includeRegex against the read whose id is matched by referenceRegex.

        @param referenceRegex: A C{str} regular expression indicating the reference
            sequence id (i.e., the sequence that others should be compared to. This must
            obviously match exactly one input sequence.
        @param window: The C{int} window size.
        @param minWindow: The C{int} minimum window size.
        @param jump: The C{int} amount to move along the genome to set up the next window.
        @param includeRegex: A regular expression indicating the sequence ids to include
            in the plot. If not specified, all sequences (other than the reference
            itself) will be compared to the reference.
        @param sort: Sort the ids of the sequences that are compared. If false, the
            order that matching ids are found in the input FASTA will be used.
        @param strict: If true, matching will be strict (ambiguous nucleotide codes
            will be considered not to match).
        @param startOffset: The C{int} zero-based genome offset to start the first
            window at.
        @param stopOffset: The C{int} zero-based genome offset to stop the last
            window at.
        @return: A C{tuple} containing the reference C{DNARead}, a C{list} of C{int}
            window start offsets, and a C{dict} keyed by C{DNARead} (for the reads whose
            ids are matched by C{includeRegex}), with values being a C{list} of C{float}
            nucleotide identity values for the windows across the genome.
        """
        self._checkReads()

        if window < 1:
            raise ValueError(f"Window size must be positive. You passed {window}.")

        if minWindow < 0:
            raise ValueError(
                f"Minimum window size cannot be negative. You passed {minWindow}."
            )

        if jump < 1:
            raise ValueError(f"Jump must be positive. You passed {jump}.")

        if startOffset < 0:
            raise ValueError(
                f"Start offset cannot be negative. You passed {startOffset}."
            )

        if stopOffset is not None:
            if stopOffset < 0:
                raise ValueError(
                    f"Stop offset cannot be negative. You passed {stopOffset}."
                )
            if stopOffset <= startOffset:
                raise ValueError(
                    f"Stop offset must be greater than start offset. You passed "
                    f"{startOffset = } and {stopOffset = }."
                )

        references = list(self._matchingReads(referenceRegex, self._reads))
        nReferences = len(references)

        # We can't use 'match' until we've dropped support for Python 3.9.
        if nReferences == 0:
            raise ValueError(
                f"No input sequence IDs match the regular expression for the reference "
                f"({referenceRegex!r})."
            )
        elif nReferences == 1:
            reference = references[0]
        else:
            ids = ", ".join(sorted(f"{read.id!r}" for read in references))
            raise ValueError(
                f"{len(references)} input sequence ids match the reference "
                f"pattern! The matching ids are: {ids}."
            )

        toCompare = list(
            set(self._matchingReads(includeRegex, self._reads)) - {reference}
        )
        nToCompare = len(toCompare)

        if nToCompare == 0:
            raise ValueError("No input sequence ids match the include pattern!")

        if sort:
            toCompare.sort(key=attrgetter("id"))

        if stopOffset is None:
            length = len(reference)
        else:
            length = min(len(reference), stopOffset)

        identity = defaultdict(list)
        starts = []
        start = startOffset

        while True:
            stop = start + window
            if stop > length:
                stop = length

            thisLength = stop - start
            if thisLength < minWindow:
                break

            starts.append(start)

            for read in toCompare:
                match = compareDNAReads(
                    reference[start:stop], read[start:stop], matchAmbiguous=not strict
                )["match"]

                identity[read].append(
                    (match["identicalMatchCount"] + match["ambiguousMatchCount"])
                    / thisLength
                )

            start += jump

        return reference, starts, identity


def addCommandLineOptions(parser):
    """
    Add options to an argument parser for any command-line tool that wants to
    work with windowed sequences.

    Note that there might be naming collisions due to the use below of simple
    option names, like --sort.
    """

    parser.add_argument(
        "--reference",
        required=True,
        help=(
            "A regular expression indicating the reference sequence id (i.e., the "
            "sequence that others should be compared to. This must obviously match "
            "exactly one input sequence."
        ),
    )

    parser.add_argument(
        "--include",
        default=".",
        help=(
            "A regular expression indicating the sequence ids to include in the plot. "
            "If not specified, all sequences (other than the reference itself) will "
            "be compared to the reference."
        ),
    )

    parser.add_argument(
        "--window",
        default=200,
        type=int,
        help="The width (number of nucleotides) of the window.",
    )

    parser.add_argument(
        "--minWindow",
        default=50,
        type=int,
        help="The smallest window size to examine.",
    )

    parser.add_argument(
        "--jump",
        default=50,
        type=int,
        help=(
            "The number of nucleotides to move across to the start of the "
            "next window."
        ),
    )

    parser.add_argument(
        "--sort",
        action="store_true",
        help=(
            "Sort the ids of the sequences that are compared. If not given, the "
            "order that matching ids are found in the input FASTA will be used."
        ),
    )

    parser.add_argument(
        "--startOffset",
        default=0,
        type=int,
        help="The zero-based genome offset to start the first window at.",
    )

    parser.add_argument(
        "--stopOffset",
        type=int,
        help="The zero-based genome offset to stop the last window at.",
    )

#!/usr/bin/env python3

import argparse
import csv
import sys

from dark.reads import (
    Read,
    ReadsInRAM,
    addFASTACommandLineOptions,
    parseFASTACommandLineOptions,
)
from dark.utils import pct


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find nucleotide sequences in an MSA and extract the matched region."
        ),
    )

    addFASTACommandLineOptions(parser)

    parser.add_argument("--prefix", help="The start of the nucleotide sequence.")

    parser.add_argument("--suffix", help="The end of the nucleotide sequence.")

    parser.add_argument(
        "--id",
        help=(
            "The id of the input sequence in which the prefix and suffix must be "
            "found. If not specified, all input sequences will be examined and "
            "those in which the prefix/suffix are found must all agree."
        ),
    )

    parser.add_argument(
        "--allowUnequalLengths",
        action="store_true",
        help=(
            "Do not exit if all sequences are not the same length. I.e., if this is "
            "not an MSA."
        ),
    )

    parser.add_argument(
        "--allowDistributedMatch",
        action="store_true",
        help=(
            "If --prefix and --suffix are both used, the default behaviour is to "
            "make sure that they are found together in at least one of the input "
            "sequences. If you don't mind if the prefix is found in one sequence "
            "and the suffix in another, use this option to allow it."
        ),
    )

    parser.add_argument(
        "--adjustIds",
        action="store_true",
        help=(
            "Add ' (position X-Y)' to the end of the output sequence ids to indicate "
            "the subsequence range of the result."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print (to standard error) a summary of total prefixes and suffixes found.",  # noqa: E501
    )

    parser.add_argument(
        "--zeroBased",
        action="store_true",
        help=(
            "If --tsv or --verbose are used, the reported prefix and suffix positions "
            "will refer to zero-based Python style [closed, open) regions."
        ),
    )

    parser.add_argument(
        "--gapCharacter", default="-", help="The 1-letter gap character."
    )

    parser.add_argument(
        "--tsv", help="The (optional) name of a TSV file to write match offsets to."
    )

    return parser.parse_args()


def summarize(
    prefix: str,
    suffix: str,
    reads: ReadsInRAM,
    overallStart: int,
    overallEnd: int,
    adjust: int,
    resultReads: ReadsInRAM,
    details: list[tuple[tuple[int, int], Read]],
    allowUnequalLengths: bool,
    gapCharacter: str,
) -> None:
    nSequences = len(reads)
    inputLengths = set(len(read) for read in reads)

    resultLengths = set(len(read) for read in resultReads)
    assert len(inputLengths) == 1 and len(resultLengths) == 1 or allowUnequalLengths
    ungappedResultLengths = set(
        len(read.sequence.replace(gapCharacter, "")) for read in resultReads
    )

    if len(inputLengths) == 1:
        length = f"{inputLengths.pop()}"
        summary = [f"Read {nSequences} sequences of length {length}."]
    else:
        lengths = ", ".join(map(str, sorted(inputLengths)))
        summary = [f"Read {nSequences} sequences (lengths {lengths}."]

    if prefix:
        summary.append(
            f"The prefix was found at position {overallStart + adjust} ({adjust}-based)."  # noqa: E501
        )
    if suffix:
        summary.append(f"The suffix ends at position {overallEnd} ({adjust}-based).")

    if len(resultLengths) == 1:
        summary.append(f"Extracted alignment region has length {resultLengths.pop()}.")
    else:
        lengths = ", ".join(map(str, sorted(resultLengths)))
        summary.append(f"Extracted aligned region sequences have lengths {lengths}.")

    if len(ungappedResultLengths) == 1:
        summary.append(
            f"Extracted ungapped sequences have length {ungappedResultLengths.pop()}."
        )
    else:
        lengths = ", ".join(map(str, sorted(ungappedResultLengths)))
        summary.append(f"Extracted ungapped sequences have lengths {lengths}.")

    prefixMatches = suffixMatches = 0
    for (start, end), _ in details:
        prefixMatches += start != -1
        suffixMatches += end != -1

    if prefix:
        summary.append(f"Prefix match count: {pct(prefixMatches, nSequences)}.")

    if suffix:
        summary.append(f"Suffix match count: {pct(suffixMatches, nSequences)}.")

    print("\n".join(summary), file=sys.stderr)


def writeTSV(
    resultReads: ReadsInRAM,
    details: list[tuple[tuple[int, int], Read]],
    tsvFile: str,
    adjust: int,
    gapCharacter: str,
) -> None:
    with open(tsvFile, "w") as fp:
        writerow = csv.writer(fp, delimiter="\t").writerow
        writerow(
            (
                "Sequence ID",
                "Prefix offset",
                "Suffix offset",
                "Match count",
                "Ungapped length",
            )
        )

        for read, ((start, end), _) in zip(resultReads, details):
            writerow(
                (
                    read.id,
                    -1 if start == -1 else start + adjust,
                    end,
                    (start != -1) + (end != -1),
                    len(read.sequence.replace(gapCharacter, "")),
                ),
            )


def main():
    args = get_args()

    prefix = args.prefix
    suffix = args.suffix
    adjust = int(not args.zeroBased)
    reads = ReadsInRAM(parseFASTACommandLineOptions(args))

    resultReads, (overallStart, overallEnd), details = reads.extractRegion(
        id_=args.id,
        prefix=prefix,
        suffix=suffix,
        ignoreGaps=True,
        allowUnequalLengths=args.allowUnequalLengths,
    )

    if overallStart != -1 and overallEnd != -1 and overallStart >= overallEnd:
        sys.exit(
            f"The prefix start position ({overallStart + adjust}) is not less "
            f"than the suffix end position ({overallEnd})"
        )

    if prefix and overallStart == -1:
        sys.exit("The prefix was not found in any sequence.")

    if suffix and overallEnd == -1:
        sys.exit("The suffix was not found in any sequence.")

    if prefix and suffix and not args.allowDistributedMatch:
        for (start, end), _ in details:
            if start != -1 and end != -1:
                # The prefix and suffix were both found in this read, so we're good.
                break
        else:
            extra = f"See {args.tsv!r}" if args.tsv else "Re-run with --tsv"
            sys.exit(
                "Simultaneous matches of the prefix and the suffix did not occur in "
                f"any one of the input sequences. {extra} for match details."
            )

    resultReads = ReadsInRAM(resultReads)

    if args.adjustIds:
        for read, ((start, end), _) in zip(resultReads, details):
            start = 0 if start == -1 else start
            end = len(read) if end == -1 else end
            read.id += f" (positions {start + adjust}-{end})"

    resultReads.save(sys.stdout)

    if args.verbose:
        summarize(
            prefix,
            suffix,
            reads,
            overallStart,
            overallEnd,
            adjust,
            resultReads,
            details,
            args.allowUnequalLengths,
            args.gapCharacter,
        )

    if args.tsv:
        writeTSV(resultReads, details, args.tsv, adjust, args.gapCharacter)


if __name__ == "__main__":
    main()

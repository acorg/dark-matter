#!/usr/bin/env python

import sys
import argparse
from math import log10
from pathlib import Path
from typing import Optional

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions, DNARead


def parseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Split a FASTA alignment into equal-sized windows.",
    )

    addFASTACommandLineOptions(parser)

    parser.add_argument(
        "--base",
        default="window-",
        help=(
            "The base output filename. If you want output to be placed in a "
            "sub-directory, specify the sub-directory in the value. Required "
            "non-existent intermediate directories (if any) will be created."
        ),
    )

    parser.add_argument("--window", type=int, default=100, help="The window size.")

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Write processing information to standard error.",
    )

    parser.add_argument(
        "--zeroBased",
        action="store_false",
        dest="oneBased",
        help=(
            "The --startOffset and offsets used in filenames will be zero-based. If "
            "not given, offsets are one-based."
        ),
    )

    parser.add_argument(
        "--minWindow",
        type=int,
        help=(
            "If the length of the final windowed sequence is less than this value, "
            "do not write a file for it."
        ),
    )

    parser.add_argument(
        "--startOffset",
        type=int,
        help=(
            "The zero-based genome offset to start windowing from (default value 0). "
            "Use --oneBased if you would prefer to use a one-based offset (in which "
            "case the default value is 1 and indices in the output files will also be "
            "one-based)."
        ),
    )

    parser.add_argument(
        "--invalidChars",
        default="N?-",
        help=(
            "Sequence characters (case insensitive) to consider invalid for the "
            "purposes of tree-making. If a windowed sequence for a sample consists "
            "entirely of these characters, no sequence will be written to the "
            "output for that window of that sample."
        ),
    )

    args = parser.parse_args()

    if args.startOffset is None:
        args.startOffset = int(args.oneBased)

    if args.oneBased:
        if args.startOffset < 1:
            sys.exit(
                "When using --oneBased, the --startOffset value must be greater than "
                "zero."
            )
    else:
        if args.startOffset < 0:
            sys.exit("The --startOffset value cannot be less than zero.")

    reads = parseFASTACommandLineOptions(args)

    return args, list(reads)


def checkReads(reads: list[DNARead]) -> int:
    """
    Check all reads have the same length and return that length.

    @param reads: A C{list} of C{DNARead} instances.
    @return: The C{int} length of the reads.
    """
    assert reads, "Empty FASTA input!"
    length = len(reads[0])
    assert length, "The first input sequence has length zero!"

    for line, read in enumerate(reads[1:], start=2):
        if len(read) != length:
            sys.exit(
                f"Sequence with id {read.id!r} on input line {line} has length "
                f"{len(read)} which is not equal to the previous input sequence "
                f"length ({length})."
            )

    return length


def writeFiles(
    reads: list[DNARead],
    window: int,
    minWindow: Optional[int],
    startOffset: int,
    oneBased: bool,
    length: int,
    base: str,
    invalidChars: str,
    verbose: bool = False,
) -> None:
    """
    Write windowed (i.e., sub-alignment) FASTA files.

    @param reads: A C{list} of C{DNARead} instances.
    @param window: The C{int} window size.
    @param minWindow: The C{int} minimum window size or C{None} if there is no minimum.
    @param startOffset: The C{int} starting offset in the genome, interpreted
        according to C{zeroBased}.
    @param oneBased: If C{True}, the passed C{startOffset} and the offsets
        used in output file names will be one-based. Else, zero-based.
    @param length: The C{int} length of all sequences in C{reads}.
    @param base: The C{str} basename of the output files.
    @param invalidChars: The C{str} (uppercase) sequence characters to ignore
        when deciding if a windowed sequence is valid.
    @param verbose: If C{True} write processing information to C{sys.stderr}.
    """
    width = int(log10(length)) + 1

    # Adjust startOffset to be zero-based.
    startOffset -= oneBased

    for start in range(startOffset, length, window):
        end = start + window
        if end > length:
            end = length

        # displayStart is the start value we print in output messages and use in
        # filenames. It is the zero-based window start offset adjusted upwards by one if
        # we have been told to use one-based values.
        displayStart = start + oneBased

        tooShort = minWindow is not None and end - start < minWindow

        filename = Path(f"{base}{displayStart:0{width}d}-{end:0{width}d}.fasta")

        if not tooShort and not filename.parent.exists():
            filename.parent.mkdir(parents=True)

        if verbose:
            print(f"Processing window {displayStart}-{end}:", file=sys.stderr)

        if tooShort:
            # Check that we are at the end of the genome.
            assert end == length
            if verbose:
                print(
                    f"  Final window (length {end - start}) shorter than the minimum window size "
                    f"({minWindow}). Skipping.",
                    file=sys.stderr,
                )
            break

        count = 0
        with open(filename, "w") as fp:
            for read in reads:
                new = read[start:end]
                sequence = new.sequence.upper()
                for char in invalidChars:
                    sequence = sequence.replace(char, "")
                if sequence:
                    count += 1
                    print(new.toString(), end="", file=fp)
                else:
                    if verbose:
                        print(
                            f"  Window {displayStart}-{end} for sequence "
                            f"{read.id!r} is all {invalidChars!r} characters. Skipping.",
                            file=sys.stderr,
                        )

        if count:
            if verbose:
                print(
                    f"  Wrote {count} window {displayStart}-{end} "
                    f"sequence{'s' if count else ''} to {str(filename)!r}.",
                    file=sys.stderr,
                )
        else:
            if verbose:
                print(
                    f"  No non-gap non-N sequences were found for window "
                    f"{displayStart}-{end}.",
                    file=sys.stderr,
                )
            filename.unlink()


def main():
    args, reads = parseArgs()
    length = checkReads(reads)
    writeFiles(
        reads,
        args.window,
        args.minWindow,
        args.startOffset,
        args.oneBased,
        length,
        args.base,
        args.invalidChars.upper(),
        args.verbose,
    )


if __name__ == "__main__":
    main()

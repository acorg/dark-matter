#!/usr/bin/env python

import sys
import argparse
from math import log10
from pathlib import Path
from typing import Optional

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions, DNARead


def getArgsAndReads():
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
        "--maxWindows",
        type=int,
        help="The maximum number of windowed subsequences to extract.",
    )

    parser.add_argument(
        "--rotate",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Rotate sequences left (if N is negative) or right (if positive) before "
            "extracting sub-sequence windows. Output filenames will be adjusted to "
            "match the rotation amount and the sequence length. For example, if the "
            "window size is 200, a --rotate value of -100 is given, and the sequences "
            "have length 2000, the first output file will have a name like "
            "window-1901-0100.fasta (or window-1900-0100.fasta if --zeroBased is "
            "used) and the final name will be window-1701-1900.fasta. This option is "
            "useful for making windows that match the location of genome features, "
            "especially in the case of a circular genome (e.g., for hepatitis B virus)."
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

    if args.maxWindows is not None and args.maxWindows < 1:
        exit("The value given to --maxWindows cannot be less than one.")

    if args.startOffset is None:
        args.startOffset = int(args.oneBased)

    if args.oneBased:
        if args.startOffset < 1:
            exit(
                "When using --oneBased, the --startOffset value must be greater than "
                "zero."
            )
    else:
        if args.startOffset < 0:
            exit("The --startOffset value cannot be less than zero.")

    reads = parseFASTACommandLineOptions(args)

    if rotate := args.rotate:
        result = [read.rotate(rotate) for read in reads]
    else:
        result = list(reads)

    return args, result


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
            exit(
                f"Sequence with id {read.id!r} on input line {line} has length "
                f"{len(read)} which is not equal to the previous input sequence "
                f"length ({length})."
            )

    return length


def writeFiles(
    reads: list[DNARead],
    window: int,
    maxWindows: int,
    minWindow: Optional[int],
    startOffset: int,
    rotate: int,
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
    @param maxWindows: The C{int} maximum number of windows to produce.
    @param minWindow: The C{int} minimum window size or C{None} if there is no minimum.
    @param startOffset: The C{int} starting offset in the genome, interpreted
        according to C{zeroBased}.
    @param rotate: The C{int} amount the reads were rotated (zero means they were not
        rotated at all). The reads have already been rotated, so this argument is only
        needed to adjust the output filenames.
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

    for windowCount, start in enumerate(range(startOffset, length, window)):

        if windowCount == maxWindows:
            break

        end = start + window
        if end > length:
            end = length

        # displayStart is the start value we print in output messages and use in
        # filenames. It is the zero-based window start offset, including any rotation
        # that has already been performed on the reads, adjusted upwards by one if we
        # have been told to use one-based values. displayEnd is similar (but does not
        # need to be adjusted by one due to Python's zero-indexing and exclusion of the
        # final index in ranges).
        displayStart = start + rotate
        displayEnd = displayStart + (end - start)

        if displayStart < 0:
            displayStart += length
        elif displayStart >= length:
            displayStart -= length

        displayStart += oneBased

        if displayEnd < 0:
            displayEnd += length
        elif displayEnd > length:
            displayEnd -= length

        tooShort = minWindow is not None and end - start < minWindow

        filename = Path(f"{base}{displayStart:0{width}d}-{displayEnd:0{width}d}.fasta")

        if not tooShort and not filename.parent.exists():
            filename.parent.mkdir(parents=True)

        if verbose:
            print(f"Processing window {displayStart}-{displayEnd}:", file=sys.stderr)

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
                            f"  Window {displayStart}-{displayEnd} for sequence "
                            f"{read.id!r} is all {invalidChars!r} characters. Skipping.",
                            file=sys.stderr,
                        )

        if count:
            if verbose:
                print(
                    f"  Wrote {count} window {displayStart}-{displayEnd} "
                    f"sequence{'s' if count else ''} to {str(filename)!r}.",
                    file=sys.stderr,
                )
        else:
            if verbose:
                print(
                    f"  No sequences with valid characters were found for window "
                    f"{displayStart}-{displayEnd}.",
                    file=sys.stderr,
                )
            filename.unlink()


def main():
    args, reads = getArgsAndReads()
    length = checkReads(reads)
    writeFiles(
        reads,
        args.window,
        args.maxWindows,
        args.minWindow,
        args.startOffset,
        args.rotate,
        args.oneBased,
        length,
        args.base,
        args.invalidChars.upper(),
        args.verbose,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python

import argparse
import multiprocessing
import sys
from dataclasses import dataclass
from typing import final

from dark.aligners import MAFFT_DEFAULT_ARGS, NEEDLE_DEFAULT_ARGS, align
from dark.fasta import FastaReads
from dark.reads import Read, Reads

MAFFT_ALGORITHMS_URL = (
    "https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html"
)


@final
class C:
    """ANSI escape codes, disabled if not a tty."""

    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    CYAN = "\033[96m"
    DIM = "\033[2m"
    BOLD = "\033[1m"
    RESET = "\033[0m"
    BG_RED = "\033[41m"
    BG_YELLOW = "\033[43m"

    @classmethod
    def disable(cls) -> None:
        for attr in [
            "RED",
            "GREEN",
            "YELLOW",
            "CYAN",
            "DIM",
            "BOLD",
            "RESET",
            "BG_RED",
            "BG_YELLOW",
        ]:
            setattr(cls, attr, "")


def pickSeq(seqs: list[Read], name: str | None, path: str) -> Read:
    """
    Select a sequence by name, or the first one if name is C{None}.

    @param seqs: A C{list} of C{Read} instances.
    @param name: A C{str} sequence id to select, or C{None} for the first.
    @param path: The C{str} file path (used in error messages).
    @return: A C{Read} instance.
    """
    if not seqs:
        print(f"Error: no sequences in {path}", file=sys.stderr)
        sys.exit(1)
    if name is None:
        return seqs[0]
    for read in seqs:
        if read.id == name:
            return read
    print(f"Error: sequence '{name}' not found in {path}", file=sys.stderr)
    print(
        f"  Available: {', '.join(r.id for r in seqs)}",
        file=sys.stderr,
    )
    sys.exit(1)


@dataclass
class Difference:
    """A single difference between two aligned sequences."""

    pos: int  # 0-based position in the alignment
    refPos: int | None  # 0-based position in seq1 (None for insertions)
    qryPos: int | None  # 0-based position in seq2 (None for deletions)
    refBase: str
    qryBase: str
    kind: str  # 'mismatch', 'deletion', 'insertion'


def findDifferences(aln1: str, aln2: str) -> list[Difference]:
    """
    Find all differences between two aligned sequences.

    @param aln1: A C{str} aligned sequence.
    @param aln2: A C{str} aligned sequence.
    @return: A C{list} of C{Difference} instances.
    """
    diffs: list[Difference] = []
    rpos, qpos = 0, 0
    for i, (a, b) in enumerate(zip(aln1, aln2)):
        if a != b:
            if a == "-":
                diffs.append(Difference(i, None, qpos, a, b, "insertion"))
            elif b == "-":
                diffs.append(Difference(i, rpos, None, a, b, "deletion"))
            else:
                diffs.append(Difference(i, rpos, qpos, a, b, "mismatch"))
        if a != "-":
            rpos += 1
        if b != "-":
            qpos += 1
    return diffs


def formatBasePair(a: str, b: str) -> tuple[str, str, str]:
    """
    Color a pair of aligned bases.

    @param a: A C{str} base from sequence 1.
    @param b: A C{str} base from sequence 2.
    @return: A 3-C{tuple} of colored C{str} values for (seq1, seq2, diff
        indicator).
    """
    if a == b:
        return C.DIM + a + C.RESET, C.DIM + b + C.RESET, " "
    elif a == "-":
        return (
            C.RED + a + C.RESET,
            C.GREEN + b + C.RESET,
            C.GREEN + "+" + C.RESET,
        )
    elif b == "-":
        return (
            C.RED + a + C.RESET,
            C.RED + b + C.RESET,
            C.RED + "-" + C.RESET,
        )
    else:
        return (
            C.YELLOW + a + C.RESET,
            C.YELLOW + b + C.RESET,
            C.YELLOW + "*" + C.RESET,
        )


def printAlignment(
    aln1: str,
    aln2: str,
    name1: str,
    name2: str,
    diffs: list[Difference],
    context: int = 10,
    width: int = 80,
) -> None:
    """
    Print alignment in blocks around differences, with context.

    @param aln1: A C{str} aligned sequence 1.
    @param aln2: A C{str} aligned sequence 2.
    @param name1: A C{str} name for sequence 1.
    @param name2: A C{str} name for sequence 2.
    @param diffs: A C{list} of C{Difference} instances.
    @param context: An C{int} number of bases of context around each
        difference.
    @param width: An C{int} alignment display width.
    """
    if not diffs:
        print(f"\n{C.GREEN}Sequences are identical ({len(aln1)} bp).{C.RESET}")
        return

    # Build regions to show: merge overlapping context windows.
    regions: list[tuple[int, int]] = []
    for d in diffs:
        start = max(0, d.pos - context)
        end = min(len(aln1), d.pos + context + 1)
        if regions and start <= regions[-1][1] + 3:
            regions[-1] = (regions[-1][0], end)
        else:
            regions.append((start, end))

    maxLabel = max(len(name1), len(name2), 3)

    for ri, (start, end) in enumerate(regions):
        if ri > 0:
            print(f"\n{C.DIM}  ...{C.RESET}\n")

        for chunkStart in range(start, end, width):
            chunkEnd = min(chunkStart + width, end)

            row1, row2, rowDiff = [], [], []
            for i in range(chunkStart, chunkEnd):
                a, b = aln1[i], aln2[i]
                fa, fb, fd = formatBasePair(a, b)
                row1.append(fa)
                row2.append(fb)
                rowDiff.append(fd)

            posLabel = f"{chunkStart + 1}"
            pad = " " * (maxLabel + 2)

            print(f"{pad}{C.DIM}{posLabel}{C.RESET}")
            print(f"  {name1:>{maxLabel}} {''.join(row1)}")
            print(f"  {'':>{maxLabel}} {''.join(rowDiff)}")
            print(f"  {name2:>{maxLabel}} {''.join(row2)}")


def printSummary(diffs: list[Difference], len1: int, len2: int, alnLen: int) -> None:
    """
    Print a summary of differences.

    @param diffs: A C{list} of C{Difference} instances.
    @param len1: The C{int} length of sequence 1.
    @param len2: The C{int} length of sequence 2.
    @param alnLen: The C{int} length of the alignment.
    """
    mismatches = sum(1 for d in diffs if d.kind == "mismatch")
    insertions = sum(1 for d in diffs if d.kind == "insertion")
    deletions = sum(1 for d in diffs if d.kind == "deletion")

    print(f"\n{C.BOLD}Summary{C.RESET}")
    print(f"  Seq 1 length:  {len1:>8,} bp")
    print(f"  Seq 2 length:  {len2:>8,} bp")
    if alnLen != len1 or alnLen != len2:
        print(f"  Alignment len: {alnLen:>8,} bp")
    print()
    print(f"  {C.YELLOW}Mismatches:  {mismatches:>6,}{C.RESET}")
    print(f"  {C.RED}Deletions:   {deletions:>6,}{C.RESET}  (in seq 1 but not seq 2)")
    print(
        f"  {C.GREEN}Insertions:  {insertions:>6,}{C.RESET}  (in seq 2 but not seq 1)"
    )
    print(f"  Total diffs:   {len(diffs):>6,}")

    matches = alnLen - len(diffs)
    identity = matches / alnLen * 100 if alnLen > 0 else 100
    print(f"  Identity:      {identity:>6.2f}%")


def printDiffTable(diffs: list[Difference], compactLimit: int = 50) -> None:
    """
    Print a tabular list of all differences.

    @param diffs: A C{list} of C{Difference} instances.
    @param compactLimit: An C{int} maximum number of differences to show
        before truncating.
    """
    if not diffs:
        return
    if len(diffs) > compactLimit:
        print(
            f"\n{C.DIM}({len(diffs)} differences — showing first "
            f"{compactLimit}){C.RESET}"
        )
        show = diffs[:compactLimit]
    else:
        print(f"\n{C.BOLD}Differences{C.RESET}")
        show = diffs

    print(f"  {'Aln':>6}  {'Seq1':>6}  {'Seq2':>6}  {'Ref':>3}  {'Qry':>3}  Type")
    print(f"  {'─' * 6}  {'─' * 6}  {'─' * 6}  {'─' * 3}  {'─' * 3}  {'─' * 10}")
    for d in show:
        rp = f"{d.refPos + 1:>6}" if d.refPos is not None else f"{'—':>6}"
        qp = f"{d.qryPos + 1:>6}" if d.qryPos is not None else f"{'—':>6}"
        kindColor = {
            "mismatch": C.YELLOW,
            "deletion": C.RED,
            "insertion": C.GREEN,
        }[d.kind]
        print(
            f"  {d.pos + 1:>6}  {rp}  {qp}"
            f"  {d.refBase:>3}  {d.qryBase:>3}"
            f"  {kindColor}{d.kind}{C.RESET}"
        )

    if len(diffs) > compactLimit:
        print(f"  {C.DIM}... and {len(diffs) - compactLimit} more{C.RESET}")


def createParser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Diff two FASTA sequences with visual highlighting.",
    )

    parser.add_argument("file1", help="First FASTA file.")
    parser.add_argument("file2", help="Second FASTA file.")

    parser.add_argument(
        "--name1",
        help="Sequence name to use from file1 (default: first).",
    )

    parser.add_argument(
        "--name2",
        help="Sequence name to use from file2 (default: first).",
    )

    parser.add_argument(
        "-c",
        "--context",
        type=int,
        default=10,
        help="Bases of context around each difference.",
    )

    parser.add_argument(
        "-w",
        "--width",
        type=int,
        default=80,
        help="Alignment display width.",
    )

    parser.add_argument(
        "--noColor",
        action="store_true",
        help="Disable colored output.",
    )

    parser.add_argument(
        "--noTable",
        action="store_true",
        help="Don't print the per-position difference table.",
    )

    parser.add_argument(
        "--maxDiffs",
        type=int,
        default=20,
        help=(
            "The maximum number of differences for which the diff table and "
            "alignment are printed. If the number of differences exceeds this "
            "value, only the summary is printed. Use -1 (the default) to "
            "always print differences."
        ),
    )

    parser.add_argument(
        "--full",
        action="store_true",
        help="Show full alignment (no context windowing).",
    )

    parser.add_argument(
        "--aligner",
        default="edlib",
        choices=("edlib", "mafft", "needle"),
        help="The alignment algorithm to use for sequences of different length.",
    )

    parser.add_argument(
        "--alignerOptions",
        help=(
            "Optional arguments to pass to the alignment algorithm. If the aligner is "
            f"'mafft', the default options are {MAFFT_DEFAULT_ARGS!r}. If 'needle', "
            f"the default is {NEEDLE_DEFAULT_ARGS!r}. Do not try to set the number of "
            "threads here; use the --threads argument instead. If you are using mafft, "
            f"see {MAFFT_ALGORITHMS_URL} for some possible option combinations."
        ),
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help=(
            "The number of threads to use when running the aligner (if the "
            "alignment algorithm can make use of multiple threads; mafft "
            "can, needle and edlib cannot)."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print alignment progress to standard error.",
    )

    return parser


def main() -> None:
    parser = createParser()

    args = parser.parse_args()

    if args.noColor or not sys.stdout.isatty():
        C.disable()

    # Read sequences.
    seqs1 = list(FastaReads(args.file1))
    seqs2 = list(FastaReads(args.file2))

    read1 = pickSeq(seqs1, args.name1, args.file1)
    read2 = pickSeq(seqs2, args.name2, args.file2)

    name1, seq1 = read1.id, read1.sequence
    name2, seq2 = read2.id, read2.sequence

    print(
        f"{C.BOLD}Comparing:{C.RESET} {name1} ({len(seq1)} bp) "
        f"vs {name2} ({len(seq2)} bp)"
    )

    # Align if different lengths, otherwise direct comparison.
    if len(seq1) == len(seq2):
        aln1, aln2 = seq1, seq2
    else:
        print(
            f"{C.DIM}Sequences differ in length, aligning with "
            f"{args.aligner}...{C.RESET}"
        )
        reads = Reads([read1, read2])
        aligned = align(
            reads,
            args.aligner,
            alignerOptions=args.alignerOptions,
            verbose=args.verbose,
            threads=args.threads,
        )
        alignedReads = list(aligned)
        aln1, aln2 = alignedReads[0].sequence, alignedReads[1].sequence

    diffs = findDifferences(aln1, aln2)
    printSummary(diffs, len(seq1), len(seq2), len(aln1))

    showDetails = args.maxDiffs == -1 or len(diffs) <= args.maxDiffs

    if showDetails and not args.noTable:
        printDiffTable(diffs)

    if showDetails:
        if args.full:
            printAlignment(
                aln1,
                aln2,
                name1,
                name2,
                diffs,
                context=len(aln1),
                width=args.width,
            )
        else:
            printAlignment(
                aln1,
                aln2,
                name1,
                name2,
                diffs,
                context=args.context,
                width=args.width,
            )


if __name__ == "__main__":
    main()

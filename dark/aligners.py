import sys
from os.path import join
from subprocess import CalledProcessError
from tempfile import mkdtemp
from shutil import rmtree
from typing import List, Optional, Tuple

import edlib

from dark.dna import AMBIGUOUS
from dark.fasta import FastaReads
from dark.process import Executor
from dark.reads import Read, Reads


EDLIB_AMBIGUOUS: Tuple[Tuple[str, ...], ...] = tuple(
    [
        tuple(sorted((nt, code)))
        for nt, ambiguities in AMBIGUOUS.items()
        for code in ambiguities
        if len(ambiguities) > 1
    ]
)


def mafft(
    reads: Reads,
    verbose: bool = False,
    options: Optional[str] = None,
    threads: Optional[int] = None,
    executor: Optional[Executor] = None,
    dryRun: bool = False,
) -> Reads:
    """
    Run a MAFFT alignment and return the sequences.

    @param reads: An iterable of multiple reads.
    @param verbose: If C{True} print progress info to sys.stderr.
    @param options: A C{str} of options to pass to mafft.
    @param executor: An C{Executor} instance. If C{None}, a new executor
        will be created and used.
    @param dryRun: If C{True}, do not execute commands, just let the executor
        log them.
    @return: A C{Reads} instance with the aligned sequences.
    """
    tempdir = mkdtemp()

    infile = join(tempdir, "input.fasta")
    out = join(tempdir, "result.fasta")

    Reads(reads).save(infile)

    if verbose:
        print("Running mafft.", file=sys.stderr)

    executor = executor or Executor(dryRun=dryRun)

    executor.execute(
        "mafft %s %s '%s' > '%s'"
        % (
            ("" if threads is None else "--thread %d" % threads),
            options or "",
            infile,
            out,
        ),
        dryRun=dryRun,
    )

    # Use 'list' in the following to force reading the FASTA from disk.
    result = Reads([] if dryRun else list(FastaReads(out)))
    rmtree(tempdir)

    return result


def needle(
    reads: List[Read],
    verbose: bool = False,
    options: Optional[str] = None,
    executor: Optional[Executor] = None,
    dryRun: bool = False,
) -> Reads:
    """
    Run a Needleman-Wunsch alignment and return the two sequences.

    @param reads: An iterable of two reads.
    @param verbose: If C{True} print progress info to sys.stderr.
    @param options: Additional options to pass to needle.
    @param executor: An C{Executor} instance. If C{None}, a new executor
        will be created and used.
    @param dryRun: If C{True}, do not execute commands, just let the executor
        log them.
    @return: A C{Reads} instance with the two aligned sequences.
    """
    tempdir = mkdtemp()

    file1 = join(tempdir, "file1.fasta")
    with open(file1, "w") as fp:
        print(reads[0].toString("fasta"), end="", file=fp)

    file2 = join(tempdir, "file2.fasta")
    with open(file2, "w") as fp:
        print(reads[1].toString("fasta"), end="", file=fp)

    out = join(tempdir, "result.fasta")

    def useStderr(e):
        return "Sequences too big. Try 'stretcher'" not in e.stderr

    executor = executor or Executor(dryRun=dryRun)

    if verbose:
        print("Running needle.", file=sys.stderr)
    try:
        executor.execute(
            "needle -asequence '%s' -bsequence '%s' %s "
            "-outfile '%s' -aformat fasta" % (file1, file2, options or "", out),
            dryRun=dryRun,
            useStderr=useStderr,
        )
    except CalledProcessError as e:
        if useStderr(e):
            raise
        else:
            if verbose:
                print(
                    "Sequences too long for needle. Falling back to "
                    "stretcher. Be patient!",
                    file=sys.stderr,
                )
            executor.execute(
                "stretcher -asequence '%s' -bsequence '%s' "
                "-auto "
                "-outfile '%s' -aformat fasta" % (file1, file2, out),
                dryRun=dryRun,
            )

    # Use 'list' in the following to force reading the FASTA from disk.
    result = Reads([] if dryRun else list(FastaReads(out)))
    rmtree(tempdir)

    return result


def removeFirstUnnecessaryGaps(
    seq1: str, seq2: str, gapSymbol: str = "-"
) -> Tuple[bool, str, str]:
    """
    Find and remove the first set of gaps in two sequences that can be removed
    without increasing the difference between the strings.

    @param seq1: A C{str} sequence string.
    @param seq2: A C{str} sequence string.
    @param gapSymbol: A C{str} 1-character symbol to use for gaps.
    @return: A 3-C{tuple} with either C{True} and two new C{str} sequences or
        C{False} and two empty strings if no removable gaps were found. The
        C{True} and C{False} first elements are just there to keep mypy happy,
        seeing as I couldn't figure out how to return (None, None) without
        getting warnings.
    """
    assert len(seq1) == len(seq2)
    for start in range(len(seq1)):
        if seq1[start] == gapSymbol:
            if seq2[start] == gapSymbol:
                raise ValueError(
                    f"Sequences {seq1!r} and {seq2!r} both have "
                    f"a gap ({gapSymbol!r}) at offset {start}!"
                )
            first, second = seq1, seq2
            break
        elif seq2[start] == gapSymbol:
            assert seq1[start] != gapSymbol
            first, second = seq2, seq1
            break
    else:
        # Did not find any gaps.
        return False, "", ""

    excessGapCount = 1
    for end in range(start + 1, len(first)):
        if first[end] == gapSymbol:
            assert second[end] != gapSymbol
            excessGapCount += 1
        elif second[end] == gapSymbol:
            excessGapCount -= 1
            if excessGapCount == 0:
                break
    else:
        # We did not get down to zero gaps. This will never happen if we
        # were originally called as a result of a call to edlibAlign
        # (below) and it only calls removeUnnecessaryGaps (also below) when
        # both original sequences have a different length from their
        # aligned versions. In that case there must be at least one gap in
        # both sequences.  But for now this is not considered an error
        # because we may not have been called by removeUnnecessaryGaps.
        return False, "", ""

    subseq1 = seq1[start : end + 1]
    subseq1noGaps = subseq1.replace(gapSymbol, "")
    subseq2 = seq2[start : end + 1]
    subseq2noGaps = subseq2.replace(gapSymbol, "")

    diffsWithGaps = sum(a != b for (a, b) in zip(subseq1, subseq2))
    diffsWithoutGaps = sum(a != b for (a, b) in zip(subseq1noGaps, subseq2noGaps))

    if diffsWithoutGaps <= diffsWithGaps:
        # The subsequences match at least as well without gaps, so replace
        # the gapped region in each sequence with its ungapped version.
        return (
            True,
            seq1[:start] + subseq1noGaps + seq1[end + 1 :],
            seq2[:start] + subseq2noGaps + seq2[end + 1 :],
        )
    else:
        return False, "", ""


def removeUnnecessaryGaps(
    seq1: str, seq2: str, gapSymbol: str = "-"
) -> Tuple[str, str]:
    """
    Find and remove all local sets of gaps in two sequences that can be removed
    without increasing the difference between the strings.

    @param seq1: A C{str} sequence string.
    @param seq2: A C{str} sequence string.
    @param gapSymbol: A C{str} 1-character symbol to use for gaps.
    @return: A 2-C{tuple} of C{str} sequences with removable gaps eliminated.
    """
    while True:
        foundGaps, new1, new2 = removeFirstUnnecessaryGaps(seq1, seq2, gapSymbol)
        if foundGaps:
            seq1, seq2 = new1, new2
        else:
            return seq1, seq2

    raise RuntimeError(
        "Fell off the end of removeFirstUnnecessaryGaps. This should be impossible!"
    )


def edlibAlign(
    reads: Reads,
    gapSymbol: str = "-",
    minimizeGaps: bool = True,
    onlyTwoSequences: bool = True,
    matchAmbiguous: bool = True,
) -> Reads:
    """
    Run an edlib alignment and return the sequences.

    @param reads: An iterable of at least two reads. The reads must not contain
        C{gapSymbol}.
    @param gapSymbol: A C{str} 1-character symbol to use for gaps.
    @param minimizeGaps: If C{True}, post-process the edlib output to remove
        unnecessary gaps.
    @param onlyTwoSequences: Ensure that C{reads} only has two reads.
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and nucleotide codes only match themselves.
    @raise ValueError: If C{onlyTwoSequences} is C{True} and there are more
       than two reads passed or if the gap symbol is already present in either
       of the reads.
    @return: A C{Reads} instance with the aligned sequences.
    """
    # Align the first two sequences.
    r1, r2, *rest = list(reads)

    # Complain if there were more than two sequences and we were told to be
    # strict.
    if onlyTwoSequences and len(rest):
        raise ValueError(f"Passed {len(rest)} unexpected extra sequences.")

    for read in r1, r2:
        if gapSymbol in read.sequence:
            raise ValueError(
                f"Sequence {read.id!r} contains one or more gap "
                f"characters {gapSymbol!r}."
            )

    alignment = edlib.getNiceAlignment(
        edlib.align(
            r1.sequence,
            r2.sequence,
            mode="NW",
            task="path",
            additionalEqualities=EDLIB_AMBIGUOUS if matchAmbiguous else None,
        ),
        r1.sequence,
        r2.sequence,
        gapSymbol=gapSymbol,
    )

    seq1, seq2 = alignment["query_aligned"], alignment["target_aligned"]

    # Try to remove unneeded gaps if requested. This is only possible if
    # the length of both sequences has changed.
    if minimizeGaps and len(seq1) != len(r1) and len(seq2) != len(r2):
        seq1, seq2 = removeUnnecessaryGaps(seq1, seq2, gapSymbol)

    return Reads([r1.__class__(r1.id, seq1), r2.__class__(r2.id, seq2)])

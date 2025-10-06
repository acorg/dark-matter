import multiprocessing
import re
import sys
from os.path import join
from shutil import rmtree
from subprocess import CalledProcessError
from tempfile import mkdtemp
from typing import Optional, Tuple

import edlib

from dark.dna import AMBIGUOUS
from dark.fasta import FastaReads
from dark.process import Executor
from dark.reads import Reads

MAFFT_DEFAULT_ARGS = "--globalpair --maxiterate 1000 --preservecase"
NEEDLE_DEFAULT_ARGS = "auto"

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
    reads: Reads,
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
        # both sequences. But for now this is not considered an error
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


# Note that this has a fixed gap char ("-"). The regex matches a non-empty sequence of
# gaps (via "-+") followed by a non-empty sequence of non-gaps (i.e., nucleotides) (via
# "[^-]+"), and collects both of the matching parts into groups named 'gaps' and 'nt'.
EDLIB_GAP_RE = re.compile("(?P<gaps>-+)(?P<nt>[^-]+)")


def _canonicalizeFirstGap(seq1: str, seq2: str, pos: int = 0) -> Tuple[int, str, str]:
    """
    Canonicalize a first gap (starting from 'pos'), if any.

    @param seq1: A C{str} sequence string.
    @param seq2: A C{str} sequence string.
    @return: A 3-C{tuple} of an int to indicate where we are up to in seq2
        (with -1 to indicate that no canonicalizations were possible)
        and two C{str} sequences with the first gap canonicalized (if possible).
    """

    for match in EDLIB_GAP_RE.finditer(seq2, pos=pos):
        gaps_start, gaps_end = match.span("gaps")
        nt_start, nt_end = match.span("nt")
        nt_len = nt_end - nt_start

        for nt_prefix_len in range(nt_len, 0, -1):
            prefix = seq2[nt_start : nt_start + nt_prefix_len]
            seq1_target = seq1[gaps_start : gaps_start + nt_prefix_len]
            if prefix == seq1_target:
                seq2 = (
                    seq2[:gaps_start]
                    + prefix
                    + "-" * (gaps_end - gaps_start)
                    + seq2[nt_start + nt_prefix_len : nt_end]
                    + seq2[nt_end:]
                )

                return gaps_start + nt_prefix_len, seq1, seq2

    # No canonicalization could be found.
    return -1, seq1, seq2


def canonicalizeAllGaps(seq1: str, seq2: str) -> Tuple[str, str]:
    """
    Canonicalize all gaps by shifting nucleotides to the right of gaps.

    @param seq1: A C{str} sequence string.
    @param seq2: A C{str} sequence string.
    @return: A 2-C{tuple} of C{str} sequences with gaps merged.
    """

    # We canonicalize by moving the nucleotides to the right (i.e., the gaps to the
    # left) because it undoes the greedy matching of edlib which matches as early as it
    # can. In many cases (with closely matching sequences where there is a gap) you
    # would get a better match (fewer gaps) if you waited. You can see an extreme
    # example of this in the final tests (e.g., testLongExample) in
    # ../../tests/test_aligners.py. That test acts on an actual output of edlib. This
    # approach to canonicalization gives you (at least in the small number of cases I
    # have looked at) the same result you get from MAFFT. MAFFT gets it "right" because
    # it has a penalty for opening a new gap, so it is disinclined to make that decision
    # if there is something more parsimonious. That also makes it much slower because it
    # has to examine many (alignment) paths. edlib just zips along doing something very
    # simple (greedy matching). Fortunately it's so basic that it's easy to improve, at
    # least in these obvious cases of just shifting the nucletotides to the right (if
    # they match at the end of the gap edlib has opened).
    #
    # The action takes place in _canonicalizeFirstGap (above), but I'm describing it
    # here where the initial and final sequence reversals happen. In between we
    # repeatedly call _canonicalizeFirstGap until it cannot find anymore gaps to
    # canonicalize.

    pos = 0
    seq1 = seq1[::-1]
    seq2 = seq2[::-1]

    while pos != -1:
        pos, seq1, seq2 = _canonicalizeFirstGap(seq1, seq2, pos)

    return seq1[::-1], seq2[::-1]


def edlibAlign(
    reads: Reads,
    gapSymbol: str = "-",
    minimizeGaps: bool = True,
    onlyTwoSequences: bool = True,
    matchAmbiguous: bool = True,
    canonicalizeGaps: bool = True,
) -> Reads:
    """
    Run an edlib alignment and return the sequences.

    @param reads: An iterable of at least two reads. The reads must not contain
        C{gapSymbol}.
    @param gapSymbol: A C{str} 1-character symbol to use for gaps.
    @param minimizeGaps: Post-process the edlib output to remove unnecessary gaps.
    @param onlyTwoSequences: Ensure that C{reads} only has two reads.
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and nucleotide codes only match themselves.
    @param canonicalizeGaps: Post-process the edlib output so that nucleotides
        that could have been at the end of a gap (instead of in the middle of
        a gap, thereby creating two gaps) are moved to the gap end. This results
        in unecessary gaps being merged.

        For example, in this alignment:

            ACTCTGACGTC
            ACT--G---TC

        the G in the second row could equivalently have been placed at the end of
        a single gap:

            ACTCTGACGTC
            ACT-----GTC

        The single gap is more akin to what MAFFT does because it is reluctant to
        open a new gap. edlib doesn't work in that way, it greedily aligns whatever
        it can as it works its way through the two sequences.

        Here's an example where gaps aren't merged but matching nucleotides are
        moved to the end of the gap. This:

            ACTCGACGTC
            ACTCG---TC

        becomes this:

            ACTCGACGTC
            ACT---CGTC

        Here is a more extreme example. This alignment with 13 gaps:

            ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC
            A-----G----AT-C-T-----G-T-----------T--C--TCT-----A-A----------A--CGAAC

        is equivalent (in terms of matched nucleotides) to this alignment with 1 gap:

            ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAAC
            --------------------------------------------------AGATCTGTTCTCTAAACGAAC

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

    if canonicalizeGaps:
        # Note that we don't pass the gap symbol here. edlib uses "-" and we can't
        # (and don't want/need to) change that.
        seq1, seq2 = canonicalizeAllGaps(seq1, seq2)

    return Reads([r1.__class__(r1.id, seq1), r2.__class__(r2.id, seq2)])


def align(
    reads: Reads,
    aligner: str,
    alignerOptions: str | None = None,
    verbose: bool = False,
    threads: int = multiprocessing.cpu_count(),
) -> Reads:
    """
    Align two reads.
    """
    if aligner == "mafft":
        options = MAFFT_DEFAULT_ARGS if alignerOptions is None else alignerOptions
        reads = mafft(reads, verbose, options=options, threads=threads)
    elif aligner == "needle":
        options = NEEDLE_DEFAULT_ARGS if alignerOptions is None else alignerOptions
        reads = needle(reads, verbose, options=options)
    else:
        assert aligner == "edlib"
        reads = edlibAlign(reads)

    return reads

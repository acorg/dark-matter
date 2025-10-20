#!/usr/bin/env python

import argparse
import csv
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from itertools import repeat
from time import time
from typing import Any, Collection

import numpy as np
import pysam

from dark.fasta import FastaReads
from dark.sam import samfile, samReferenceLengths
from dark.utils import baseCountsToStr, entropy2, gmt, pct

BASES = "ACGT"
TRANSITIONS = "AG GA CT TC".split()
TRANSVERSIONS = "AT TA AC CA GT TG GC CG".split()


def makeParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Print SAM/BAM per-site nucleotide composition and coverage statistics. "
            "Note that no output is produced for sites with no read coverage."
        ),
    )

    parser.add_argument(
        "samfile",
        help="The SAM/BAM file to read.",
    )

    parser.add_argument(
        "--referenceId",
        metavar="ID",
        help=(
            "A reference sequence id. If not given and the SAM/BAM file has just one "
            "reference, that one will be used. Otherwise, a reference id must be given."
        ),
    )

    parser.add_argument(
        "--fasta",
        help=(
            "Optionally give a FASTA file in which the reference can be found. "
            "If you do not provide a reference sequence, mutations will be counted "
            "as differences from the most-common nucleotide at each genome site."
        ),
    )

    parser.add_argument(
        "--diffsFrom",
        choices=("reference", "commonest"),
        default="reference",
        help=(
            "The base from which to assess differences. If 'reference', differences "
            "from the reference base (at each genome site) will be treated as "
            "mutations. If 'commonest', differences from the most common base will be "
            "treated as mutations (with ties broken in favour of the reference base, "
            "if a reference is given). It is important that you understand the "
            "difference between these two and choose the right option."
        ),
    )

    parser.add_argument(
        "--window",
        "-n",
        type=int,
        default=250,
        help=(
            "The number of genome sites to have each subprocess handle. Large values "
            "will have a major performance impact when --maxDepth is also high (e.g., "
            ">10,000). This seems to be due to an unknown pysam/htslib slowdown (based "
            "on samtools v1.21, as of 2025-01-31)."
        ),
    )

    parser.add_argument(
        "--workers",
        type=int,
        help=(
            "The number of subprocesses to run in parallel. If not specified, the "
            "process pool will typically use all available cores, but see "
            "https://docs.python.org/3/library/concurrent.futures.html"
            "#processpoolexecutor for details."
        ),
    )

    parser.add_argument(
        "--noOffsets",
        action="store_false",
        dest="printOffsets",
        help="Do not print per-site details of base counts.",
    )

    parser.add_argument(
        "--noStats",
        action="store_false",
        dest="printStats",
        help="Do not print final average, standard deviation etc. statistics.",
    )

    parser.add_argument(
        "--statsFile",
        help=(
            "The file to print overall statistics to. If not specified, standard "
            "output will be used."
        ),
    )

    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help=(
            "Print progress information (to standard error). 0 = don't print progress "
            "information. 1 = print information after each genome region is processed. "
            "2 = also print information during the processing of each region."
        ),
    )

    parser.add_argument(
        "--minDepth",
        type=int,
        help=(
            "The minimum read depth to report on. Sites with a lower depth will "
            "not be reported."
        ),
    )

    parser.add_argument(
        "--maxDepth",
        type=int,
        default=1e6,
        help="The maximum read depth to request from pysam for each site.",
    )

    parser.add_argument(
        "--minSiteMutationRate",
        type=float,
        help=(
            "The minimum site mutation rate to report on. Sites with a lower mutation "
            "rate will not be reported. The mutation rate is the fraction of bases at "
            "the site that are not the reference or commonest base (depending on "
            "the --diffsFrom option)."
        ),
    )

    parser.add_argument(
        "--maxSiteMutationRate",
        type=float,
        help=(
            "The maximum site mutation rate to report on. Sites with a higher mutation "
            "rate will not be reported. The mutation rate is the fraction of bases at "
            "the site that are not the reference or commonest base (depending on "
            "the --diffsFrom option)."
        ),
    )

    parser.add_argument(
        "--minSiteMutationCount",
        type=int,
        default=0,
        help=(
            "The minimum number of mutant nucleotides that must be present at a site "
            "in order for it to be reported. E.g., if this were set to 4, there would "
            "need to be at least four nucleotides (of any mixture of bases) that were "
            "not the consensus or reference (depending on the --diffsFrom setting)."
        ),
    )

    parser.add_argument(
        "--minEntropy",
        type=float,
        default=0.0,
        help=(
            "The minimum site entropy to report on. Sites with lower entropy "
            "will not be reported."
        ),
    )

    parser.add_argument(
        "--maxEntropy",
        type=float,
        default=2.0,
        help=(
            "The maximum site entropy to report on. Sites with higher entropy "
            "will not be reported."
        ),
    )

    parser.add_argument(
        "--noFilter",
        action="store_false",
        dest="filter",
        help="This behaviour is now the default and this option is ignored.",
    )

    parser.add_argument(
        "--tsv",
        action="store_true",
        help=(
            "Produce TSV output (use --sep , if you want CSV). Note that overall "
            "SAM/BAM file statistics will be printed at the bottom of the output "
            "unless you also specify --noStats or --statsFile. Columns will hopefully "
            "be self-explanatory, but a few clarifying comments will probably be "
            "helpful. 1) The 'Base' column indicates the nucleotide base that was used "
            "to count mutations. The base used will depend on the value of --diffsFrom. "  # noqa: E501
            "See the help output for that option. 2) The 'Reference' column will have "
            "an empty value if no reference sequence was provided. 3) The 'Site' column "  # noqa: E501
            "gives 1-based indexes into the reference genome."
        ),
    )

    parser.add_argument(
        "--sep",
        default="\t",
        help="The separator to use when --tsv is used. Ignored otherwise.",
    )

    parser.add_argument(
        "--minBases",
        type=int,
        choices=(1, 2, 3, 4),
        default=1,
        help=(
            "The minimum number of different bases that must occur at a site for "
            "the site to be reported. E.g., if this were set to 3, sites are "
            "only counted if at most one nucleotide is not present. The default "
            "value of 1 will do no filtering because all (covered) sites will have "
            "at least one nucleotide present. Higher values of this parameter can "
            "be used to add a diversity filter, by only counting sites at which "
            "there is greater nucleotide diversity."
        ),
    )

    parser.add_argument(
        "--maxIdenticalReadIdsPerSite",
        type=int,
        default=2,
        help=(
            "The number of reads with identical IDs that are permitted at each site. "
            "This should be 2 for paired-end reads (Illumina sequencing uses the "
            "identical read ID for both members of the pair) because the read pair "
            "may overlap at some sites. If you have unpaired reads, you could set this "
            "to 1. A factor that needs to be considered is whether the mapping "
            "algorithm was told to only produce the best match for each read (or read "
            "pair). If so, setting this to 1 or 2 (as just described) makes sense "
            "because it will catch inadvertent errors. But if reads are allowed to "
            "match multiple times, there could be many matches of the same read in "
            "the same reference section. In that case, you can use a value of -1 to "
            "disable the checking."
        ),
    )

    return parser


def getReferenceId(
    filename: str, references: Collection[str], referenceId: str | None
) -> str:
    """
    Get the reference id.  This is either the reference id from the command line
    (if given and if present in the SAM/BAM file) or the single reference from the
    SAM/BAM file (if it contains just one reference).

    @param filename: The SAM/BAM filename.
    @param references: The names of all references in the SAM/BAM file.
    @param referenceId: The id of the wanted reference, if provided on the command line.
    @return: The reference id.
    """
    if referenceId:
        if referenceId not in references:
            raise ValueError(
                f"Reference id {referenceId} not found in {filename!r} "
                f"Known references are {', '.join(sorted(references))}."
            )
    else:
        if len(references) == 1:
            referenceId = list(references)[0]
        else:
            raise ValueError(
                f"Multiple references found in {filename!r}: "
                f"{', '.join(sorted(references))}. Use --referenceId "
                "to specify the one you want."
            )

    return referenceId


def getReferenceSeq(filename: str, referenceId: str) -> str:
    """
    Look for the reference id in the FASTA file and return its sequence if found.
    If no reference is found, raise a ValueError.

    @param filename: The FASTA filename.
    @param referenceId: The ID of the wanted reference.
    @return: The reference nucleotide sequence.
    """
    reads = list(FastaReads(filename))

    # First look for matches with the full id+description sequence.
    for read in reads:
        if read.id == referenceId:
            return read.sequence

    # Then look for matches with the first (space separated) part of the full
    # id+description sequence. We do this second so as to allow the full id+description
    # to be attempted first and only if that doesn't work do we allow for a
    # less-complete match.
    for read in reads:
        if read.id.split()[0] == referenceId:
            return read.sequence

    raise ValueError(
        f"Reference id {referenceId!r} was not found in reference FASTA "
        f"file {filename!r}."
    )


def analyzeColumn(
    column: pysam.PileupColumn,
    start: int,
    window: int,
    minDepth: int | None,
    minSiteMutationRate: float | None,
    maxSiteMutationRate: float | None,
    minSiteMutationCount: int,
    minEntropy: float,
    maxEntropy: float,
    minBases: int,
    maxIdenticalReadIdsPerSite: int,
    verbosity: int,
    referenceSeq: str | None,
    diffsFrom: str,
) -> (
    tuple[
        int,
        set[str],
        float,
        float,
        dict[str, int],
        int,
        str | None,
        str,
        dict[str, int],
    ]
    | None
):
    """
    Do an analysis of the reads in a column. Return a tuple of results
    unless the column should be skipped, in which case return None.

    @param start: The start offset of the genome window we process.
    @param window: The size of the genome section we should process.
    @param miDepth: The minimum read depth to report on. Sites with lower depth will
        be ignored.
    @param minSiteMutationRate: The minimum mutation rate (i.e., fraction) required at
        sites in order that they are reported. Sites with lower mutation fractions will
        be ignored.
    @param maxSiteMutationRate: The maximum mutation rate (i.e., fraction) allowed at
        sites in order that they are reported. Sites with higher mutation fractions will
        be ignored.
    @param minSiteMutationCount: The minimum number of mutant nucleotides that must be
        present at a site in order for it to be reported. E.g., if this were set to 4,
        there would need to be at least four nucleotides (of any mixture of bases) that
        were not the consensus or reference (depending on the --diffsFrom setting).
    @param minEntropy: The minimum entropy to allow.
    @param maxEntropy: The maximum entropy to allow.
    @param minBases: The minimum number of different bases that must occur at a site for
        the site to be reported. See the parser argument above for more details (or run
        this script with --help).
    @param maxIdenticalReadIdsPerSite: The number of reads with identical IDs that are
        permitted at each site. Use -1 to turn off the checking. See the parser argument
        above for more details (or run this script with --help).
    @param verbosity: The verbosity level.
    @param referenceSeq: The genome of the reference (if known).
    @param diffsFrom: The base from which to assess differences. If 'reference',
        differences from the reference base (at each genome site) will be treated as
        mutations. If 'commonest', differences from the most common base will be
        treated as mutations (with ties broken in favour of the reference base,
        if a reference is given).
    @return: A tuple containing the following variables (types are above in the
        typehint), whose names are hopefully self-explanatory:

            read ids
            mutation rate
            entropy
            mutation pair counts
            reference offset
            reference base
            from base
            bases
    """
    refOffset = column.reference_pos
    assert start <= refOffset < start + window

    # readIds will have read ids as keys, with values that are lists containing:
    #
    #   the number of times the read has been seen in this column (site).
    #   the current highest quality for the base for that read at this column (site).
    #   the base with the best quality.
    #
    # This information is collected because both reads of a read pair can
    # overlap the same site.  We want to choose the base that has the highest
    # quality and to detect if there is somehow more than a single read pair
    # (with the same read id) matching a site (this could happen if a read
    # mapping algorithm maps the same read pair more than once and the
    # regions of the two matches overlaps).
    readIds: dict[str, list[Any]] = {}

    for read in column.pileups:
        # assert isinstance(read, pysam.PileupRead)
        if read.query_position is None:
            assert read.is_del or read.is_refskip
        else:
            # assert isinstance(read.alignment, pysam.AlignedSegment)
            readId = read.alignment.query_name
            assert readId
            assert read.alignment.query_sequence is not None  # Typing check.
            base = read.alignment.query_sequence[read.query_position]

            if base not in BASES:
                # In practice this should never happen because the FASTQ we analyze
                # always only has ACGT bases. Allowing non-ACGT bases through here
                # results in a KeyError below when we count base frequencies using a
                # dict.fromkeys(BASES, 0). So issue a warning and continue.
                print(
                    f"WARNING: Site {refOffset + 1} is matched by a read "
                    f"(id {readId!r}) containing a non-ACGT base ({base}). "
                    "This is very unusual, as FASTQ data should normally only "
                    "contain sequences composed of unambiguous ACGT nucleotide calls.",
                    file=sys.stderr,
                )
                continue

            assert read.alignment.query_qualities is not None  # Typing check.
            quality: str = read.alignment.query_qualities[read.query_position]
            try:
                count, previousQuality, previousBase = readIds[readId]
            except KeyError:
                # We've not seen this read at this site before. Remember the
                # details so we can compare against them later (in the 'else'
                # just below) if the read pair also overlaps this site.
                readIds[readId] = [1, quality, base]
            else:
                count += 1
                if (
                    maxIdenticalReadIdsPerSite > 0
                    and count > maxIdenticalReadIdsPerSite
                ):
                    raise ValueError(
                        f"Read id {readId} occurs (at least) {count} times for site "
                        f"{refOffset + 1}. The maximum allowed repeats of the same read "  # noqa: E501
                        f"ID at a site is {maxIdenticalReadIdsPerSite}."
                    )
                # This is either the read pair or another match of the same
                # read to this site.
                readIds[readId][0] = count

                if base != previousBase:
                    if quality > previousQuality:
                        readIds[readId][1] = quality
                        readIds[readId][2] = base
                        if verbosity:
                            print(
                                f"Site {refOffset + 1}: replaced base "
                                f"({previousBase!r} with {base!r}) due to improved "
                                "quality in second read.",
                                file=sys.stderr,
                            )
                    elif quality == previousQuality:
                        # We should probably choose the consensus (or
                        # reference, if known) base here, if they are among
                        # the options. But we don't know what the consensus
                        # base is, yet. This may not be a problem as a later
                        # read (with the same ID) might have a base with a
                        # higher quality.
                        if verbosity > 1:
                            print(
                                f"Site {refOffset + 1}: has multiple bases "
                                f"({base!r} and {previousBase!r}) of the same quality. "
                                f"Ignoring {base!r}. Note that we may encounter "
                                "another read for this site with a higher quality, in "
                                "which case these equally-lower-quality base calls "
                                "will be ignored in its favour.",
                                file=sys.stderr,
                            )

    nReads = sum(readInfo[0] for readInfo in readIds.values())

    if nReads == 0 or (minDepth is not None and nReads < minDepth):
        # nReads can still be zero at this point because all reads
        # matching this site have a read.query_position of None (due
        # to having deletion or reference skips).
        if verbosity > 1:
            print(
                f"Site {refOffset + 1}: skipped due to too few reads ({nReads:,}).",
                file=sys.stderr,
            )
        return

    bases = dict.fromkeys(BASES, 0)
    for _, _, base in readIds.values():
        bases[base] += 1

    # Set an initial nonsense entropy value to keep type checking from
    # warning that the variable may be unbound.
    entropy = -1.0

    nt = []
    for base, count in bases.items():
        nt.extend([base] * count)
        entropy = entropy2(nt)

    if entropy < minEntropy:
        if verbosity > 1:
            print(
                f"Site {refOffset + 1}: skipped due to low entropy ({entropy:.5f}).",
                file=sys.stderr,
            )
        return

    if entropy > maxEntropy:
        if verbosity > 1:
            print(
                f"Site {refOffset + 1}: skipped due to high entropy ({entropy:.5f}).",
                file=sys.stderr,
            )
        return

    refBase = referenceSeq[refOffset] if referenceSeq else None

    if diffsFrom == "reference":
        assert refBase is not None
        fromBase = refBase
    else:
        # Find the base with the highest count. This is the 'from' base
        # (i.e., the originally existing base that was mutated from something
        # to something else). Ties are broken in favour of the reference base
        # (if there is no reference, the second value in the compared tuples
        # will always be False (because refBase will be None)). If there are
        # two equally frequent bases, neither of which is the reference (if
        # known), the tie is broken in favour of the alphabetically first
        # base.
        fromBase = max(
            (count, int(base == refBase), base) for (base, count) in bases.items()
        )[2]

    mutationRate = 1.0 - bases[fromBase] / nReads

    if minSiteMutationRate is not None and mutationRate < minSiteMutationRate:
        if verbosity > 1:
            print(
                f"Site {refOffset + 1}: skipped due to low mutation rate "
                f"({mutationRate:.5f}).",
                file=sys.stderr,
            )
        return

    if maxSiteMutationRate is not None and mutationRate > maxSiteMutationRate:
        if verbosity > 1:
            print(
                f"Site {refOffset + 1}: skipped due to high mutation rate "
                f"({mutationRate:.5f}).",
                file=sys.stderr,
            )
        return

    mutationCount = 0
    mutationPairCounts: dict[str, int] = defaultdict(int)
    for base in BASES:
        if base != fromBase and (count := bases[base]):
            mutationPairCounts[fromBase + base] += count
            mutationCount += count

    if minSiteMutationCount > 0 and mutationCount < minSiteMutationCount:
        if verbosity > 1:
            print(
                f"Site {refOffset + 1}: skipped due to having insufficient "
                f"nucleotide diversity ({mutationCount} different nucleotides "
                "are present, but at least {minSiteMutationCount} are required).",
                file=sys.stderr,
            )
        return

    # The number of bases at the site is the 'from' base, plus the bases it
    # has mutated to. Hence the + 1 in the following test.
    if len(mutationPairCounts) + 1 < minBases:
        if verbosity > 1:
            print(
                f"Site {refOffset + 1}: skipped due to having insufficient "
                f"nucleotide diversity ({len(mutationPairCounts) + 1} different "
                f"nucleotides are present, but at least {minBases} are required).",
                file=sys.stderr,
            )
        return

    return (
        nReads,
        set(readIds),
        mutationRate,
        entropy,
        mutationPairCounts,
        refOffset,
        refBase,
        fromBase,
        bases,
    )


def makeRow(
    refOffset: int,
    refBase: str | None,
    nReads: int,
    bases: dict[str, int],
    fromBase: str,
    mutationRate: float,
    mutationPairCounts: dict[str, int],
    entropy: float,
) -> list[int | str | float | None]:
    """
    Make a row of data (for a site) to be added to the rows returned by
    C{processWindowColumns}.

    @param refOffset: The 0-based offset into the reference.
    @param refBase: The reference base (if the reference is known).
    @param nReads: The number of reads for the column being examined.
    @param bases: The count of each nucleotide base at this site.
    @param fromBase: The base that is considered canonical at this site (i.e.,
        the one that mutations are measured from).
    @param mutationRate: The fraction of bases at the site that are mutations.
    @param mutationPairCounts: The count of each mutation pair at this site, where keys
        are the concatenated 'from' and 'to' bases involved in the mutation (e.g.,
        mutationPairCounts["CT"] holds the number of C -> T mutations).
    @param entropy: The entropy of the bases at this genome offset.
    """
    mutationCount = sum(mutationPairCounts.values())

    row = (
        [refOffset + 1, refBase, nReads]
        # Don't use bases.values() in the following, because we need bases in
        # ACGT order.
        + [bases[base] for base in BASES]
        + [
            fromBase,
            len(mutationPairCounts) + 1,
            mutationCount,
            round(mutationRate, 6),
            round(entropy, 6),
        ]
    )

    totalTransitions = 0
    for pair in TRANSITIONS:
        count = mutationPairCounts.get(pair, 0)
        row.append(count)
        totalTransitions += count

    totalTransversions = 0
    for pair in TRANSVERSIONS:
        count = mutationPairCounts.get(pair, 0)
        row.append(count)
        totalTransversions += count

    row.extend((totalTransitions, totalTransversions))

    return row


def processWindowColumns(
    filename: str,
    start: int,
    window: int,
    minDepth: int | None,
    maxDepth: int,
    minSiteMutationRate: float | None,
    maxSiteMutationRate: float | None,
    minSiteMutationCount: int,
    minEntropy: float,
    maxEntropy: float,
    minBases: int,
    maxIdenticalReadIdsPerSite: int,
    collectOffsets: bool,
    verbosity: int,
    referenceId: str,
    referenceSeq: str | None,
    returnList: bool,
    startTime: float,
    diffsFrom: str,
) -> dict[str, Any]:
    """
    Process all columns in a subsection (window) of the genome.

    @param filename: The SAM/BAM file to read.
    @param start: The start offset of the genome window we process.
    @param window: The size of the genome section we should process.
    @param miDepth: The minimum read depth to report on. Sites with lower depth will
        be ignored.
    @param maxDepth: The maximum read depth to consider.
    @param minSiteMutationRate: The minimum mutation rate (i.e., fraction) required at
        sites in order that they are reported. Sites with lower mutation fractions will
        be ignored.
    @param maxSiteMutationRate: The maximum mutation rate (i.e., fraction) allowed at
        sites in order that they are reported. Sites with higher mutation fractions will
        be ignored.
    @param minSiteMutationCount: The minimum number of mutant nucleotides that must be
        present at a site in order for it to be reported. E.g., if this were set to 4,
        there would need to be at least four nucleotides (of any mixture of bases) that
        were not the consensus or reference (depending on the --diffsFrom setting).
    @param minBases: The minimum number of different bases that must occur at a site for
        the site to be reported. See the parser argument above for more details (or run
        this script with --help).
    @param maxIdenticalReadIdsPerSite: The number of reads with identical IDs that are
        permitted at each site. Use -1 to turn off the checking. See the parser argument
        above for more details (or run this script with --help).
    @param collectOffsets: Whether to collect genome offset information.
    @param verbosity: The verbosity level.
    @param referenceId: The name of the reference sequence.
    @param referenceSeq: The genome of the reference (if known).
    @param returnList: If True, the rows we return will be printed (by our caller) as
        CSV, so we return lists. Else, we return a string for each row of output.
    @param startTime: The time overall processing was launched.
    @param diffsFrom: The base from which to assess differences. If 'reference',
        differences from the reference base (at each genome site) will be treated as
        mutations. If 'commonest', differences from the most common base will be
        treated as mutations (with ties broken in favour of the reference base,
        if a reference is given).
    """
    startTime = time()
    allReadIds = set()
    rows = []
    depths = []
    cumulativeReadCount = mutationCount = 0
    baseCounts = dict.fromkeys(BASES, 0)
    totalMutationPairCounts = defaultdict(int)

    with samfile(filename) as sam:
        for column in sam.pileup(
            start=start,
            end=start + window,
            reference=referenceId,
            truncate=True,
            max_depth=maxDepth,
        ):
            assert isinstance(column, pysam.PileupColumn)  # Typing check.
            columnStartTime = time()

            columnInfo = analyzeColumn(
                column,
                start,
                window,
                minDepth,
                minSiteMutationRate,
                maxSiteMutationRate,
                minSiteMutationCount,
                minEntropy,
                maxEntropy,
                minBases,
                maxIdenticalReadIdsPerSite,
                verbosity,
                referenceSeq,
                diffsFrom,
            )

            if columnInfo is None:
                # This column should be skipped. It didn't have enough diversity or
                # reads, or failed some other check.
                continue

            (
                nReads,
                columnReadIds,
                mutationRate,
                entropy,
                mutationPairCounts,
                refOffset,
                refBase,
                fromBase,
                columnBaseCounts,
            ) = columnInfo

            if diffsFrom == "commonest":
                assert columnBaseCounts[fromBase] >= max(columnBaseCounts.values())

            depths.append(nReads)
            allReadIds.update(columnReadIds)
            mutationCount += sum(mutationPairCounts.values())

            for base, count in columnBaseCounts.items():
                baseCounts[base] += count

            for pair, count in mutationPairCounts.items():
                totalMutationPairCounts[pair] += count

            if collectOffsets:
                if returnList:
                    rows.append(
                        makeRow(
                            refOffset,
                            refBase,
                            nReads,
                            columnBaseCounts,
                            fromBase,
                            mutationRate,
                            mutationPairCounts,
                            entropy,
                        )
                    )
                else:
                    out = f"{refOffset + 1} {nReads} {baseCountsToStr(baseCounts)}"
                    if referenceSeq:
                        out += f" Ref:{refBase}"
                    rows.append(out)

            if verbosity > 1:
                cumulativeReadCount += nReads
                stopTime = time()
                elapsed = stopTime - columnStartTime
                totalElapsed = stopTime - startTime

                print(
                    f"Site {refOffset + 1}: {nReads:,} reads processed "
                    f"at {nReads / elapsed / 1e6:.2f} M read/s. "
                    f"Total bases: {sum(baseCounts.values()):,} "
                    f"Overall: {cumulativeReadCount / totalElapsed / 1e6:.2f} M read/s",
                    file=sys.stderr,
                )

    return {
        "depths": depths,
        "baseCounts": baseCounts,
        "mutationPairCounts": totalMutationPairCounts,
        "readIds": allReadIds,
        "rows": rows,
        "startSite": start,
        "startTime": startTime,
        "stopTime": time(),
        "mutationCount": mutationCount,
    }


def printStats(
    statsFile,
    samfile,
    minDepth,
    maxDepth,
    diffsFrom,
    minBases,
    maxIdenticalReadIdsPerSite,
    referenceSeq,
    referenceId,
    referenceLen,
    startGMT,
    elapsed,
    readIds,
    depths,
    minSiteMutationRate,
    maxSiteMutationRate,
    minSiteMutationCount,
    minEntropy,
    maxEntropy,
    totalBaseCounts,
    mutationPairCounts,
) -> None:
    """
    Print a summary of the run. See the docstrings elsewhere if you want to know what
    all these parameters mean :-)
    """
    with open(statsFile, "w") if statsFile else sys.stdout as out, redirect_stdout(out):
        print("Inputs:")
        print("  SAM/BAM file:", samfile)
        print("  Reference:")
        print("    ID:", referenceId)
        print("    Length:", referenceLen)
        print("    Sequence given:", "yes" if referenceSeq else "no")
        print("Settings:")
        print("  Max identical read IDs per site:", maxIdenticalReadIdsPerSite)
        print("  Depth:")
        print("    Min depth required:", minDepth)
        print("    Max depth considered:", maxDepth)
        print("  Mutation (by site):")
        print("    Counts are relative to:", diffsFrom)
        print("    Min mutation rate:", minSiteMutationRate)
        print("    Max mutation rate:", maxSiteMutationRate)
        print(f"    Min entropy: {minEntropy:.5f}")
        print(f"    Max entropy: {maxEntropy:.5f}")
        print(f"    Min non-{diffsFrom} nucleotides:", minSiteMutationCount)
        print("    Min different nucleotides:", minBases)
        print("Result:")
        print("  Timing:")
        print("    Started at:", startGMT)
        print("    Stopped at:", gmt())
        print(f"    Elapsed: {elapsed} seconds")
        print(
            f"  Number of reads covering sites of interest: {sum(depths)} "
            f"({len(readIds)} unique read ids)"
        )
        print("  Reference genome sites included:", pct(len(depths), referenceLen))
        print("  Coverage depth:")
        print("    Min:", min(depths, default=0))
        print("    Max:", max(depths, default=0))

        if depths:
            print(f"    Mean: {np.mean(depths):.3f}")
            print(f"    Median: {np.median(depths):.3f}")
            print(f"    SD: {np.std(depths):.3f}")

        print("  Bases:")
        totalBases = 0
        for base, count in totalBaseCounts.items():
            print(f"    {base}: {count}")
            totalBases += count
        print(f"    TOTAL: {totalBases}")

        mutationCount = sum(mutationPairCounts.values())
        overallMutationRate = mutationCount / totalBases if totalBases else 0.0

        print("  Mutations:")
        print(f"    Overall rate: {overallMutationRate:.6f}")

        totalTransitions = totalTransversions = 0
        print("    Sorted counts (* = transversion):")
        countWidth = None

        def key(pair):
            return mutationPairCounts[pair], pair

        for pair in sorted(mutationPairCounts, key=key, reverse=True):
            count = mutationPairCounts[pair]
            if countWidth is None:
                countWidth = len(str(count))
            if pair in TRANSITIONS:
                totalTransitions += count
                what = ""
            else:
                assert pair in TRANSVERSIONS
                totalTransversions += count
                what = " *"
            print(f"      {pair[0]} -> {pair[1]}: {count:{countWidth}d}{what}")

        print("    Transitions:", totalTransitions)
        print("    Transversions:", totalTransversions)
        print("    TOTAL:", totalTransitions + totalTransversions)


def main() -> None:
    args = makeParser().parse_args()

    if not (args.printOffsets or args.printStats):
        sys.exit("You used both --noOffsets and --noStats, so there is no output!")

    if args.maxIdenticalReadIdsPerSite == 0:
        sys.exit("maxIdenticalReadIdsPerSite cannot be zero.")

    referenceLens = samReferenceLengths(args.samfile)
    referenceId = getReferenceId(args.samfile, set(referenceLens), args.referenceId)
    referenceLen = referenceLens[referenceId]
    referenceSeq = getReferenceSeq(args.fasta, referenceId) if args.fasta else None

    if args.diffsFrom == "reference" and referenceSeq is None:
        sys.exit(
            "If you use --diffsFrom 'reference' (the default), you must use --fasta to "
            "give a FASTA file containing the reference sequence. Else use --diffsFrom "
            "'commonest' to calculate mutations based on the most-common nucleotide "
            "at each genome site."
        )

    printTSV = csv.writer(sys.stdout, delimiter=args.sep).writerow if args.tsv else None
    printRow = printTSV or print

    if printTSV:
        # Print the TSV header.
        printTSV(
            ("Site", "Reference", "Depth")
            + tuple(BASES)
            + ("Base", "nBases", "Mutations", "Variability", "Entropy")
            + tuple(f"{pair[0]}->{pair[1]}" for pair in TRANSITIONS + TRANSVERSIONS)
            + ("Transitions", "Transversions"),
        )

    depths = []
    readIds = set()
    mutationCount = 0
    # I don't use a collections.Counter in the following because they are very slow.
    baseCounts = defaultdict(int)
    mutationPairCounts = defaultdict(int)
    start = time()
    startGMT = gmt()
    totalRunTime = maxRunTime = 0.0
    startSites = list(range(0, referenceLen, args.window))

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for result in executor.map(
            processWindowColumns,
            repeat(args.samfile),
            startSites,
            repeat(args.window),
            repeat(args.minDepth),
            repeat(args.maxDepth),
            repeat(args.minSiteMutationRate),
            repeat(args.maxSiteMutationRate),
            repeat(args.minSiteMutationCount),
            repeat(args.minEntropy),
            repeat(args.maxEntropy),
            repeat(args.minBases),
            repeat(args.maxIdenticalReadIdsPerSite),
            repeat(args.printOffsets),
            repeat(args.verbosity),
            repeat(referenceId),
            repeat(referenceSeq),
            repeat(bool(printTSV)),
            repeat(start),
            repeat(args.diffsFrom),
        ):
            waited = result["startTime"] - start
            runTime = result["stopTime"] - result["startTime"]
            totalRunTime += runTime
            if runTime > maxRunTime:
                maxRunTime = runTime

            if args.verbosity > 1:
                from_ = result["startSite"]
                to = min(from_ + args.window, referenceLen)
                print(
                    f"Genome region {from_ + 1}-{to} processed in {runTime:.2f} "
                    f"seconds (waited {waited:.2f} seconds to start).",
                    file=sys.stderr,
                )

            if args.printOffsets:
                for row in result["rows"]:
                    printRow(row)

            readIds.update(result["readIds"])
            depths.extend(result["depths"])
            for base, count in result["baseCounts"].items():
                baseCounts[base] += count
            for pair, count in result["mutationPairCounts"].items():
                mutationPairCounts[pair] += count
            mutationCount += result["mutationCount"]

    stop = time()
    elapsed = int(stop - start)

    if args.verbosity > 0:
        try:
            speedup = totalRunTime / (stop - start)
        except ZeroDivisionError:
            speedup = 0.0

        print(
            f"Finished in {elapsed} seconds. Total serial (i.e., non parallel) "
            f"processing would have been {int(totalRunTime)} seconds. The slowest "
            f"process took {maxRunTime:.2f} seconds. Speed up = {speedup:.2f}",
            file=sys.stderr,
        )

    if args.printStats:
        printStats(
            args.statsFile,
            args.samfile,
            args.minDepth,
            args.maxDepth,
            args.diffsFrom,
            args.minBases,
            args.maxIdenticalReadIdsPerSite,
            referenceSeq,
            referenceId,
            referenceLen,
            startGMT,
            elapsed,
            readIds,
            depths,
            args.minSiteMutationRate,
            args.maxSiteMutationRate,
            args.minSiteMutationCount,
            args.minEntropy,
            args.maxEntropy,
            baseCounts,
            mutationPairCounts,
        )


if __name__ == "__main__":
    main()

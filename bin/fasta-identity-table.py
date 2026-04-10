#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict
from operator import itemgetter
from typing import Iterable

from dark.aligners import edlibAlign, mafft, needle
from dark.dna import compareDNAReads
from dark.filter import (
    addFASTAEditingCommandLineOptions,
    addFASTAFilteringCommandLineOptions,
    parseFASTAEditingCommandLineOptions,
    parseFASTAFilteringCommandLineOptions,
)
from dark.reads import (
    Read,
    addFASTACommandLineOptions,
    getNoCoverageCounts,
    parseFASTACommandLineOptions,
)

MAFFT_DEFAULT_ARGS = "--globalpair --maxiterate 1000 --preservecase"
MAFFT_ALGORITHMS_URL = (
    "https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html"
)
NEEDLE_DEFAULT_ARGS = "auto"

# What to put in an empty HTML table cell.
EMPTY_CELL = "&nbsp;"


def align(reads: Iterable[Read], args: argparse.Namespace) -> list[Read]:
    """
    Align a pair of reads.

    @param reads: An iterable of two C{DNARead} instances.
    @param args: An argparse C{Namespace} instance with command-line options.
    @return: A C{list} of two aligned C{DNARead} instances.
    """
    if args.aligner == "mafft":
        # Be careful in examining args.alignerOptions because we want the
        # user to be able to pass an empty string (so check against None
        # before deciding to use the default.)
        options = (
            MAFFT_DEFAULT_ARGS if args.alignerOptions is None else args.alignerOptions
        )
        return mafft(reads, args.verbose, options=options, threads=args.threads)
    elif args.aligner == "needle":
        # Be careful in examining args.alignerOptions because we want the
        # user to be able to pass an empty string (so check against None
        # before deciding to use the default.)
        options = (
            NEEDLE_DEFAULT_ARGS if args.alignerOptions is None else args.alignerOptions
        )
        return needle(reads, args.verbose, options=options)
    else:
        assert args.aligner == "edlib"
        return edlibAlign(reads)


def thresholdToCssName(threshold: float) -> str:
    """
    Turn a floating point threshold into a string that can be used as a CSS
    class name.

    param threshold: The C{float} threshold.
    return: A C{str} CSS class name.
    """
    # '.' is illegal in a CSS class name.
    return "threshold-" + str(threshold).replace(".", "_")


def thresholdForIdentity(identity: float, colors: list[tuple[float, str]]) -> float:
    """
    Get the best identity threshold for a specific identity value.

    @param identity: A C{float} nucleotide identity.
    @param colors: A C{list} of (threshold, color) tuples, where threshold is a
        C{float} and color is a C{str} to be used as a cell background. This
        is as returned by C{parseColors}.
    @return: The first C{float} threshold that the given identity is at least
        as big as.
    """
    for threshold, _ in colors:
        if identity >= threshold:
            return threshold
    raise ValueError("This should never happen! Last threshold is not 0.0?")


def parseColors(colors: list[str], defaultColor: str) -> list[tuple[float, str]]:
    """
    Parse command line color information.

    @param colors: A C{list} of space separated "value color" strings, such as
        ["0.9 red", "0.75 rgb(23, 190, 207)", "0.1 #CF3CF3"].
    @param defaultColor: The C{str} color to use for cells that do not reach
        the identity fraction threshold of any color in C{colors}.
    @return: A C{list} of (threshold, color) tuples, where threshold is a
        C{float} (from C{colors}) and color is a C{str} (from C{colors}). The
        list will be sorted by decreasing threshold values.
    """
    result = []
    if colors:
        for colorInfo in colors:
            fields = colorInfo.split(None, 1)
            if len(fields) == 2:
                threshold, color = fields
                try:
                    threshold = float(threshold)
                except ValueError:
                    sys.exit(
                        "--color arguments must be given as space-separated "
                        'pairs of "value color" where the value is a numeric '
                        f"identity threshold. Your value {threshold!r} is not numeric."
                    )

                if 0.0 > threshold > 1.0:
                    sys.exit(
                        "--color arguments must be given as space-separated "
                        'pairs of "value color" where the value is a '
                        "numeric identity threshold from 0.0 to 1.0. Your "
                        f"value {threshold!r} is outside that range."
                    )

                result.append((threshold, color))
            else:
                sys.exit(
                    "--color arguments must be given as space-separated pairs of "
                    f'"value color". You have given {colorInfo!r}, which does not '
                    "contain a space."
                )

    result.sort(key=itemgetter(0), reverse=True)

    if not result or result[-1][0] > 0.0:
        result.append((0.0, defaultColor))

    return result


def getGapCounts(reads: Iterable[Read], gapChars: str) -> dict[str, int]:
    """
    Get the number of gap characters in all reads.

    @param reads: A C{Reads} instance.
    @param gapChars: A C{str} of sequence characters considered to be gaps.
    @return: A C{dict} keyed by read id, with C{int} gap counts.
    """
    return {read.id: sum(c in set(gapChars) for c in read.sequence) for read in reads}


def explanation(args: argparse.Namespace) -> str:
    """
    Make an explanation of the output HTML table.

    @param args: An argparse C{Namespace} instance with command-line options.
    @return: A C{str} of HTML.
    """
    result = [
        """
<h1>Sequence versus sequence identity table</h1>

<p>

The table cells below show the nucleotide identity fraction for the sequences
(<span class="best">like this</span> for the best value in each row). The
identity fraction numerator is the sum of the number of identical
    """
    ]

    if args.matchAmbiguous:
        result.append(
            "nucleotides plus the number of ambiguously matching nucleotides."
        )
    else:
        result.append("nucleotides.")

    result.append(
        """The denominator
is the length of the sequence <em>for the row</em>. Sequence gaps
are not included when calculating their lengths.

</p>
    """
    )

    if (
        args.showLengths
        or args.showGaps
        or args.showNs
        or args.matchAmbiguous
        or args.showMatchDetails
    ):
        result.append(
            """
<p>

Key to abbreviations:
  <ul>
    """
        )

        if args.showLengths:
            result.append("<li>L: sequence <strong>L</strong>ength.</li>")

        if args.showGaps:
            result.append("<li>G: number of <strong>G</strong>aps in sequence.</li>")

        if args.showNoCoverage:
            result.append(
                "<li>C: number of no <strong>C</strong>overage "
                "characters in sequence.</li>"
            )

        if args.showNs:
            result.append(
                "<li>N: number of fully-ambiguous "
                "<strong>N</strong> characters in sequence.</li>"
            )

        if args.showMatchDetails:
            result.append(
                "<li>IM: <strong>I</strong>dentical nucleotide "
                "<strong>M</strong>atches.</li>"
            )

        if args.matchAmbiguous:
            result.append(
                "<li>AM: <strong>A</strong>mbiguous nucleotide "
                "<strong>M</strong>atches.</li>"
            )

        if args.showGaps:
            result.append(
                """
    <li>GG: <strong>G</strong>ap/<strong>G</strong>ap matches
                (both sequences have gaps).</li>
    <li>G?: <strong>G</strong>ap/Non-gap mismatches (one sequence has a gap).</li>
"""
            )

        if args.showNoCoverage:
            result.append(
                """
    <li>CC: No <strong>c</strong>overage/No <strong>c</strong>overage
                (both sequences have no coverage).</li>
    <li>C?: No <strong>c</strong>overage (one sequence has no coverage).</li>
"""
            )

    result.append(
        """
    <li>NE: <strong>N</strong>on-<strong>e</strong>qual nucleotide mismatches.</li>
  </ul>
</p>
"""
    )

    return "\n".join(result)


def dataCell(
    rowId: str,
    colId: str,
    symmetric: bool,
    rowReadIndices: dict[str, int],
    colReadIndices: dict[str, int],
    upperOnly: bool,
) -> bool:
    """
    Should a cell have a value computed for it, or is it empty?

    @param rowId: The row read id.
    @param colId: The column read id.
    @param symmetric: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param rowReadIndices: A C{dict} mapping read ids to 0-based row indices.
    @param colReadIndices: A C{dict} mapping read ids to 0-based column indices.
    @param upperOnly: If C{True}, only compute values for the upper diagonal.
    @return: C{True} if the cell should be filled in, else C{False}.
    """
    if rowId == colId and symmetric:
        return False

    return rowReadIndices[rowId] < colReadIndices[colId] if upperOnly else True


def collectData(
    rowReads: dict[str, Read],
    colReads: dict[str, Read],
    symmetric: bool,
    args: argparse.Namespace,
) -> tuple[dict[str, dict[str, int | list[int]]], dict[str, int], dict[str, int]]:
    """
    Get pairwise matching statistics for two sets of reads.

    @param rowReads: An C{dict} of C{str} read ids whose values are
        C{Read} instances. These will be the rows of the table.
    @param colReads: An C{dict} of C{str} read ids whose values are
        C{Read} instances. These will be the columns of the table.
    @param symmetric: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param args: An argparse C{Namespace} instance with command-line options.
    @return: A 3-tuple containing:
        - A C{dict} whose keys are read (row) ids, with values being a C{dict} whose
          keys are read (column) ids with values being the DNA match (a C{dict},
          returned by C{compareDNAReads}) between the two reads.
        - A C{dict} mapping C{str} read ids to 0-based row indices.
        - A C{dict} mapping C{str} read ids to 0-based column indices.
    """
    rowReadIndices = {readId: n for n, readId in enumerate(rowReads)}
    colReadIndices = {readId: n for n, readId in enumerate(colReads)}
    comparisons = sum(
        dataCell(
            rowId, colId, symmetric, rowReadIndices, colReadIndices, args.upperOnly
        )
        for rowId in rowReads
        for colId in colReads
    )

    result = defaultdict(dict)
    count = 0

    for rowId, rowRead in rowReads.items():
        for colId, colRead in colReads.items():
            if dataCell(
                rowId, colId, symmetric, rowReadIndices, colReadIndices, args.upperOnly
            ):
                count += 1
                if args.align:
                    r1, r2 = align([rowRead, colRead], args)
                else:
                    r1, r2 = rowRead, colRead
                if args.verbose:
                    print(
                        f"Comparing {count}/{comparisons} {rowId!r} and {colId!r}.",
                        file=sys.stderr,
                    )
                match = compareDNAReads(
                    r1,
                    r2,
                    matchAmbiguous=args.matchAmbiguous,
                    gapChars=args.gapChars,
                    noCoverageChars=args.noCoverageChars,
                )
                if not args.matchAmbiguous:
                    assert match["match"]["ambiguousMatchCount"] == 0
                # Record the lengths, since these may have changed due to making the
                # alignment.
                match["read1"]["length"] = len(r1)
                match["read2"]["length"] = len(r2)
                result[rowId][colId] = result[colId][rowId] = match

    return result, rowReadIndices, colReadIndices


def computeIdentity(
    rowRead: Read,
    colRead: Read,
    stats: dict[str, dict[str, int | list[int]]],
    matchAmbiguous: bool,
    digits: int,
) -> float:
    """
    Compute nucleotide identity for two reads (as a fraction of the lowest number of
    relevant nucleotides in either read).

    @param rowRead: A C{Read} instance.
    @param colRead: A C{Read} instance.
    @param stats: A C{dict} as returned by C{compareDNAReads}.
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are possibly
        correct as actually being correct. Otherwise, we are strict and insist that
        only non-ambiguous nucleotides can contribute to the matching nucleotide count.
    @param digits: The C{int} number of digits to round values to.

    """

    # Note that the strict identity may be higher or lower than the ambiguous identity
    # even though an ambiguous match sounds like it should always be more lenient and
    # therefore result in a better match. While it is true that the raw number of
    # matched characters will always increase when strict=False, the fraction of matched
    # characters may go down.
    #
    # The strict value can be higher because rowRead might have many ambiguous
    # characters but very few of them may match colRead. In that case the overall
    # fraction of matching characters will be pulled down from the strict fraction when
    # the ambiguous are included.
    #
    # Similarly, rowRead may have many ambiguous characters, all of which are matched by
    # colRead and this can pull the overall identity higher than the strict identity.

    match = stats["match"]
    denominator = min(stats["read1"]["length"], stats["read2"]["length"]) - (
        match["gapGapMismatchCount"]
        + match["noCoverageCount"]
        + match["noCoverageNoCoverageCount"]
    )
    numerator = stats["match"]["identicalMatchCount"]
    if matchAmbiguous:
        numerator += stats["match"]["ambiguousMatchCount"]

    result = numerator / denominator
    assert result <= 1.0, f"{numerator} / {denominator} = {result}.\n{stats}"

    return round(result, digits)


def textTable(
    tableData,
    rowReads: dict[str, Read],
    colReads: dict[str, Read],
    rowReadIndices: dict[str, int],
    colReadIndices: dict[str, int],
    symmetric: bool,
    args: argparse.Namespace,
) -> None:
    """
    Print (to standard output) a text table showing inter-sequence distances.

    @param tableData: A C{defaultdict(dict)} keyed by read ids, whose values
        are the dictionaries returned by compareDNAReads.
    @param rowReads: An C{dict} of C{str} read ids whose values are
        C{Read} instances. These will be the rows of the table.
    @param colReads: An C{dict} of C{str} read ids whose values are
        C{Read} instances. These will be the columns of the table.
    @param rowReadIndices: A C{dict} mapping read ids to 0-based row indices.
    @param colReadIndices: A C{dict} mapping read ids to 0-based column indices.
    @param symmetric: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param args: An argparse C{Namespace} instance with command-line options.
    """
    titles = ["ID"]
    if args.numberedColumns:
        titles.extend(str(i + 1) for i in range(len(colReads)))

        if args.upperOnly:
            titles.pop(1)
            titles[-1] = list(colReads)[-1]
    else:
        titles.extend(colReads)

    print("\t".join(titles))

    for rowCount, (rowId, rowRead) in enumerate(rowReads.items(), start=1):
        if rowCount == len(rowReads) and args.concise:
            # We don't print the last row for a concise table, because it will be the
            # final column with a value against all rows, so the final row can be
            # omitted.
            pass
        prefix = f"{rowCount}: " if args.numberedColumns else ""
        print(f"{prefix}{rowId}", end="")
        for colId, colRead in colReads.items():
            if colReadIndices[colId] == 0 and args.concise:
                # The whole first column is skipped if we're making a concise table
                # because the values for that column are in the first row.
                continue
            if dataCell(
                rowId, colId, symmetric, rowReadIndices, colReadIndices, args.upperOnly
            ):
                identity = computeIdentity(
                    rowRead,
                    colRead,
                    tableData[rowId][colId],
                    args.matchAmbiguous,
                    args.digits,
                )

                if args.addZeroes:
                    print(f"\t{identity:.{args.digits}f}", end="")
                else:
                    print(f"\t{identity}", end="")
            else:
                print("\t", end="")
        print()


def htmlTable(
    tableData,
    rowReads: dict[str, Read],
    colReads: dict[str, Read],
    rowReadIndices: dict[str, int],
    colReadIndices: dict[str, int],
    symmetric: bool,
    colors: list[tuple[float, str]],
    args: argparse.Namespace,
) -> str:
    """
    Make an HTML table showing inter-sequence distances.

    @param tableData: A C{defaultdict(dict)} keyed by read ids, whose values
        are the dictionaries returned by compareDNAReads.
    @param rowReads: An C{dict} of C{str} read ids whose values are
        C{Read} instances. These will be the rows of the table.
    @param colReads: An C{dict} of C{str} read ids whose values are
        C{Read} instances. These will be the columns of the table.
    @param rowReadIndices: A C{dict} mapping read ids to 0-based row indices.
    @param colReadIndices: A C{dict} mapping read ids to 0-based column indices.
    @param symmetric: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param colors: A C{list} of (threshold, color) tuples, where threshold is a
        C{float} and color is a C{str} to be used as a cell background. This
        is as returned by C{parseColors}.
    @param args: An argparse C{Namespace} instance with command-line options.
    @return: An HTML C{str} showing inter-sequence distances.
    """
    gaps1 = getGapCounts(rowReads.values(), args.gapChars)
    gaps2 = getGapCounts(colReads.values(), args.gapChars)
    noCoverage1 = getNoCoverageCounts(rowReads.values(), args.noCoverageChars)
    noCoverage2 = getNoCoverageCounts(colReads.values(), args.noCoverageChars)
    result = []
    append = result.append
    extend = result.extend

    def writeHeader():
        # The header row of the table.
        extend(
            [
                "    <tr>",
                f'    <td class="empty">{EMPTY_CELL}</td>',
            ]
        )

        for count, colRead in enumerate(colReads.values(), start=1):
            if count == 1 and args.concise:
                # The first column will be skipped, so skip it in the header too.
                continue
            append(
                '    <td class="title"><span class="name">%s</span>'
                % (count if args.numberedColumns else colRead.id)
            )
            if not symmetric:
                if args.showLengths:
                    append("    <br>L:%d" % len(colRead))
                if args.showGaps:
                    append("    <br>G:%d" % gaps2[colRead.id])
                if args.showNoCoverage:
                    append("    <br>C:%d" % noCoverage2[colRead.id])
                if args.showNs:
                    append("    <br>N:%d" % colRead.sequence.count("N"))
            append("    </td>")
        append("    </tr>")

    if args.div:
        append("<div>")
    else:
        extend(
            [
                "<!DOCTYPE HTML>",
                "<html>",
                "<head>",
                '<meta charset="UTF-8">',
                "</head>",
                "<body>",
            ]
        )

    extend(
        [
            "<style>",
            f"""
        table {{
            border-collapse: collapse;
        }}
        table, td {{
            border: 1px solid #ccc;
        }}
        tr:hover {{
            background-color: #f2f2f2;
        }}
        td {{
            vertical-align: top;
            font-size: 14px;
        }}
        td.empty {{
            min-width: {args.emptyCellMinSide};
            min-height: {args.emptyCellMinSide};
        }}
        span.name {{
            font-weight: bold;
        }}
        span.best {{
            font-weight: bold;
        }}
        td.nt-identity {{
            text-align: right;
        }}
                """,
        ]
    )

    # Add color style information for the identity thresholds.
    for threshold, color in colors:
        append(".%s { background-color: %s; }" % (thresholdToCssName(threshold), color))

    append("</style>")

    if not args.div:
        append(explanation(args))

    extend(
        [
            '<div style="overflow-x:auto;">',
            "<table>",
            "  <tbody>",
        ]
    )

    # Pre-process to find the best identities in each sample row.
    bestIdentityForId = {}
    identities = defaultdict(dict)

    for rowId, rowRead in rowReads.items():
        # Look for best identity for the sample.
        bestIdentity = -1.0
        for colId, colRead in colReads.items():
            if dataCell(
                rowId, colId, symmetric, rowReadIndices, colReadIndices, args.upperOnly
            ):
                identity = computeIdentity(
                    rowRead,
                    colRead,
                    tableData[rowId][colId],
                    args.matchAmbiguous,
                    args.digits,
                )
                identities[rowId][colId] = identity
                if identity > bestIdentity:
                    bestIdentity = identity

        bestIdentityForId[rowId] = bestIdentity

    writeHeader()

    # The main body of the table.
    for rowCount, (rowId, rowRead) in enumerate(rowReads.items(), start=1):
        if rowCount == len(rowReads) and args.concise:
            # We don't print the last row for a concise table, because it will be the
            # final column with a value against all rows, so the final row can be
            # omitted.
            continue

        extend(
            [
                "    <tr>",
                '      <td class="title"><span class="name">%s%s</span>'
                % (f"{rowCount}: " if args.numberedColumns else "", rowId),
            ]
        )
        if args.showLengths:
            append("<br/>L:%d" % len(rowRead))
        if args.showGaps:
            append("<br/>G:%d" % gaps1[rowRead.id])
        if args.showNoCoverage:
            append("<br/>C:%d" % noCoverage1[rowRead.id])
        if args.showNs:
            append("<br/>N:%d" % rowRead.sequence.count("N"))
        append("</td>")

        for colId, colRead in colReads.items():
            if colReadIndices[colId] == 0 and args.concise:
                # The whole first column is skipped if we're making a concise table
                # because the values for that column are in the first row.
                continue

            if not dataCell(
                rowId, colId, symmetric, rowReadIndices, colReadIndices, args.upperOnly
            ):
                append(f'<td class="empty">{EMPTY_CELL}</td>')
                continue

            identity = identities[rowId][colId]

            append(
                '      <td class="nt-identity %s">'
                % thresholdToCssName(thresholdForIdentity(identity, colors))
            )

            # The maximum percent identity.
            if args.highlightBest and identity == bestIdentityForId[rowId]:
                scoreStyle = ' class="best"'
            else:
                scoreStyle = ""

            if args.addZeroes:
                append(f"<span{scoreStyle}>{identity:.{args.digits}f}</span>")
            else:
                append(f"<span{scoreStyle}>{identity}</span>")

            if args.showMatchDetails:
                match = tableData[rowId][colId]["match"]
                append("<br/>IM:%d" % match["identicalMatchCount"])

                if args.matchAmbiguous:
                    append("<br/>AM:%d" % match["ambiguousMatchCount"])

                if args.showGaps:
                    append(
                        "<br/>GG:%d<br/>G?:%d"
                        % (match["gapGapMismatchCount"], match["gapMismatchCount"])
                    )

                if args.showNoCoverage:
                    append(
                        "<br/>CC:%d<br/>C?:%d"
                        % (match["noCoverageCount"], match["noCoverageNoCoverageCount"])
                    )

                append("<br/>NE:%d" % match["nonGapMismatchCount"])
            append("      </td>")
        append("    </tr>")

    if args.footer:
        writeHeader()

    extend(
        [
            "  </tbody>",
            "</table>",
            "</div>",
        ]
    )

    if args.div:
        append("</div>")
    else:
        extend(
            [
                "</body>",
                "</html>",
            ]
        )

    return "\n".join(result)


def main():
    parser = argparse.ArgumentParser(
        description="Print a FASTA sequence identity table."
    )

    parser.add_argument(
        "--text",
        action="store_true",
        help="If specified, just print a simple text table.",
    )

    parser.add_argument(
        "--strict",
        action="store_false",
        dest="matchAmbiguous",
        help="If given, do not allow ambiguous nucleotide symbols to match.",
    )

    parser.add_argument(
        "--numberedColumns",
        action="store_true",
        help=(
            "Use a sequence (row) number as the header of each column "
            "instead of the sequence id."
        ),
    )

    parser.add_argument(
        "--omitMatchDetails",
        action="store_false",
        dest="showMatchDetails",
        help="Do not show match details.",
    )

    parser.add_argument(
        "--digits",
        default=3,
        type=int,
        help="The number of digits to round identities to.",
    )

    parser.add_argument(
        "--addZeroes",
        action="store_true",
        help=(
            "If given, make sure all identities have the same number of "
            "digits shown by adding zeroes when needed."
        ),
    )

    parser.add_argument(
        "--showLengths",
        action="store_true",
        help="If given, show the lengths of sequences.",
    )

    parser.add_argument(
        "--showGaps",
        action="store_true",
        help="If given, show the number of gaps in sequences.",
    )

    parser.add_argument(
        "--showNoCoverage",
        action="store_true",
        help=("If given, show the number of no-coverage characters in sequences."),
    )

    parser.add_argument(
        "--showNs",
        action="store_true",
        help=(
            "If given, show the number of fully ambiguous N characters in sequences."
        ),
    )

    parser.add_argument(
        "--footer",
        action="store_true",
        help="If given, also show sequence ids at the bottom of the table.",
    )

    parser.add_argument(
        "--div",
        action="store_true",
        help=(
            "If given, print an HTML <div> fragment only, not a full HTML "
            "document (ignored if --text is used)."
        ),
    )

    parser.add_argument(
        "--fastaFile2",
        metavar="FILENAME",
        help=(
            "The name of a second FASTA input file. If no second FASTA "
            "file is given, sequences in the first FASTA file will be "
            "compared with each other."
        ),
    )

    parser.add_argument(
        "--gapChars",
        default="-",
        metavar="CHARS",
        help=(
            "The sequence characters that should be considered to be gaps. "
            "These characters will be ignored in computing sequence lengths "
            "and identity fractions."
        ),
    )

    parser.add_argument(
        "--noCoverageChars",
        metavar="CHARS",
        help=(
            "The sequence characters that indicate lack of coverage. "
            "These characters will be ignored in identity fractions."
        ),
    )

    parser.add_argument(
        "--defaultColor",
        default="white",
        help=(
            "The (background) color for cells. This will be used for all "
            "cells that do not otherwise have a color due to use of "
            "--color. This option is ignored if --text is given."
        ),
    )

    parser.add_argument(
        "--color",
        action="append",
        help=(
            "Specify cell background coloring. This option must be given as "
            'a space-separated "value color" pair. The value is an identity '
            "fraction in [0..1] and the color is any color "
            "specification that can be given to CSS. This argument can be "
            'repeated. E.g., --color "0.9 red" --color "0.75 rgb(23, 190, '
            '207)" --color "0.1 #CF3CF3". Cells will be colored using the '
            "color of the highest identity fraction they satisfy. The "
            "default is to color all cells with the --defaultColor color. "
            "This option is ignored if --text is given."
        ),
    )

    parser.add_argument(
        "--align", action="store_true", help="Do pairwise alignment of all sequences."
    )

    parser.add_argument(
        "--upperOnly",
        action="store_true",
        help=(
            "Only show the upper diagonal (only valid when a single FASTA "
            "input file is given, resulting in a square all-against-all "
            "identity table)."
        ),
    )

    parser.add_argument(
        "--concise",
        action="store_true",
        help=(
            "When only one file is given (and the table is therefore square with "
            "symmetric entries), Make the table more concise by omitting the first "
            "column and last row. These can be omitted because the first column has "
            "the sequence that is in the first row (where it is matched against "
            "everything else) and the last row is the same sequence as the final "
            "column (where it is matched against everything else). When these two "
            "are omitted, the table is slightly smaller and arguably slightly harder "
            "to understand."
        ),
    )

    parser.add_argument(
        "--sort", action="store_true", help="Sort the input sequences by id."
    )

    parser.add_argument(
        "--aligner",
        default="edlib",
        choices=("edlib", "mafft", "needle"),
        help="The alignment algorithm to use.",
    )

    parser.add_argument(
        "--emptyCellMinSide",
        default="70px",
        help="The CSS minimum width and height for empty HTML table cells.",
    )

    parser.add_argument(
        "--alignerOptions",
        help=(
            "Optional arguments to pass to the alignment algorithm. If the "
            'aligner is mafft, the default options are %r. If needle, "%s". '
            "Do not try to set the number of threads here - use the "
            "--threads argument instead. If you are using mafft, see %s "
            "for some possible option combinations."
            % (MAFFT_DEFAULT_ARGS, NEEDLE_DEFAULT_ARGS, MAFFT_ALGORITHMS_URL)
        ),
    )

    parser.add_argument(
        "--verbose", action="store_true", help="Print progress to standard error."
    )

    parser.add_argument(
        "--highlightBest",
        action="store_true",
        help="Highlight the best identity in each row.",
    )

    addFASTACommandLineOptions(parser)
    addFASTAFilteringCommandLineOptions(parser)
    addFASTAEditingCommandLineOptions(parser)

    args = parser.parse_args()

    if args.concise and args.numberedColumns:
        sys.exit(
            "If you use --concise, you cannot also use --numberedColumns. Othewise, "
            "the number at the top of the final column cannot be translated into a "
            "sequence ID since the final row is not printed in a concise table."
        )

    colors = parseColors(args.color, args.defaultColor)
    # Sanity check - the last threshold must be zero.
    assert colors[-1][0] == 0.0

    reads = parseFASTAEditingCommandLineOptions(
        args,
        parseFASTAFilteringCommandLineOptions(args, parseFASTACommandLineOptions(args)),
    )

    # Collect the reads into a dict, keeping the insertion order, unless we
    # are told to sort.
    rowReads = {}
    for read in sorted(reads) if args.sort else reads:
        rowReads[read.id] = read

    if args.fastaFile2:
        if args.upperOnly:
            sys.exit(
                "The --upperOnly option is not supported when two FASTA input files "
                "are given."
            )
        if args.concise:
            sys.exit(
                "The --concise option is not supported when two FASTA input files "
                "are given."
            )

        symmetric = False

        _tmpReads = {}
        # The next line is a total hack, to trick parseFASTACommandLineOptions
        # into reading a second FASTA file.
        args.fastaFile = args.fastaFile2
        for read in parseFASTAFilteringCommandLineOptions(
            args, parseFASTACommandLineOptions(args)
        ):
            _tmpReads[read.id] = read
        colReads = {}
        for readId in sorted(_tmpReads) if args.sort else _tmpReads:
            colReads[readId] = _tmpReads[readId]
        del _tmpReads
    else:
        symmetric = True
        colReads = rowReads

    tableData, rowReadIndices, colReadIndices = collectData(
        rowReads,
        colReads,
        symmetric,
        args,
    )

    if args.text:
        textTable(
            tableData,
            rowReads,
            colReads,
            rowReadIndices,
            colReadIndices,
            symmetric,
            args,
        )
    else:
        print(
            htmlTable(
                tableData,
                rowReads,
                colReads,
                rowReadIndices,
                colReadIndices,
                symmetric,
                colors,
                args,
            )
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python

import sys
import argparse
from operator import itemgetter
from itertools import cycle, chain
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objs as go
from plotly.io import write_image, write_html
from collections import Counter, defaultdict

from dark.windowedIdentity import WindowedIdentity, addCommandLineOptions
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions, DNARead


def makeParser() -> argparse.ArgumentParser:
    """
    Make an argument parser.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Make a plot showing windowed nucleotide identity across an alignment."
        ),
    )

    addFASTACommandLineOptions(parser)
    addCommandLineOptions(parser)

    parser.add_argument(
        "--html",
        help="The name of an HTML file to write a plotly image to.",
    )

    parser.add_argument(
        "--out",
        help=(
            "The name of the file to write the image to. The file suffix "
            "will determine the output format. The image will be created with "
            "plotly unless --matplotlib is specified."
        ),
    )

    parser.add_argument(
        "--title",
        help="The plot title.",
    )

    parser.add_argument(
        "--legendTitle",
        default="Sample",
        help="The legend title.",
    )

    parser.add_argument(
        "--dryRun",
        "-n",
        action="store_true",
        help="Just print what would be compared without actually doing it.",
    )

    parser.add_argument(
        "--bestOnly",
        action="store_true",
        help="Only show the best matching sequences for each window.",
    )

    parser.add_argument(
        "--matplotlib",
        action="store_true",
        help="Use matplotlib to make a static image. Otherwise, plotly will be used.",
    )

    parser.add_argument(
        "--y01",
        action="store_true",
        help="Set the y-axis limits to (0.0, 1.0). Only applicable for matplotlib plots.",
    )

    parser.add_argument(
        "--ambiguous",
        dest="strict",
        action="store_false",
        help=(
            "Count compatible ambiguous nucleotide codes as matching when computing "
            "identity."
        ),
    )

    parser.add_argument(
        "--markerSize",
        type=int,
        help="The size of the markers in plot lines.",
    )

    return parser


def plotMatplotlib(
    centers: list[float],
    identity: dict[DNARead, list[float]],
    title: str,
    args: argparse.Namespace,
) -> None:
    """
    Make a matplotlib image.

    @param centers: The x-axis centers of the genome windows.
    @param identity: A mapping from sequences (that have been matched against the
        reference) to a list of float identities (one for each window).
    @param title: The title of the plot.
    @param args: The command-line arguments.
    """
    figure, ax = plt.subplots(1, 1, figsize=(10, 8))

    for read, identities in identity.items():
        ax.plot(centers, identities, label=f"{read.id}")

    if args.y01:
        ax.set_ylim((0.0, 1.0))

    ax.grid(axis="x", color="0.95")
    ax.legend(title=args.legendTitle)
    figure.suptitle(title)

    figure.savefig(args.out)


def getBestOnly(
    centers: list[float], identity: dict[DNARead, list[float]], args: argparse.Namespace
) -> list[go.Scatter]:
    """
    Find the best-matching sequence(s) for each window and return scatter plots for them.

    @param centers: The x-axis centers of the genome windows.
    @param identity: A mapping from sequences (that have been matched against the
        reference) to a list of float identities (one for each window).
    @param title: The title of the plot.
    @param args: The command-line arguments.
    """
    tolerance = 0.99
    delta = 0.005
    half = args.jump / 2.0
    inLegend = set()

    mode, marker = (
        ("lines+markers", dict(size=args.markerSize))
        if args.markerSize  # i.e., not 0 and not None.
        else ("lines", None)
    )

    # Keep track of how many times each read has a window that is amongst the
    # best. We'll use that to update the legend to display the count.
    bestCounts: Counter[str] = Counter()

    # scatters will hold a list of go.Scatter instances, keyed by readId.
    scatters = defaultdict(list)

    colors: dict[str, str] = {}
    colorIterator = cycle(
        chain(px.colors.qualitative.Dark2, px.colors.qualitative.Pastel2)
    )

    for jumpIndex in range(len(centers)):
        windowIdentities = []
        # It's less efficient to loop over .items() repeatedly like this, but I think it
        # makes for simpler code since we process each window in turn.
        for read, identities in identity.items():
            windowIdentities.append((identities[jumpIndex], read.id))
        # Sort first by ascending read id.
        windowIdentities.sort(key=itemgetter(1))
        # Then by descending identity.
        windowIdentities.sort(key=itemgetter(0), reverse=True)

        bestMatches = [windowIdentities[0]]
        bestIdentity = bestMatches[0][0]

        # Collect all sequences that are within the tolerance of the best identity.
        for thisIdentity, readId in windowIdentities[1:]:
            if thisIdentity >= tolerance * bestIdentity:
                bestMatches.append((thisIdentity, readId))
            else:
                break

        center = centers[jumpIndex]
        x = [center - half, center + half]

        # Collect the best matches.
        lastIdentity = 2.0  # Must be > 1.0 to start things off.
        for thisIdentity, readId in bestMatches:
            bestCounts[readId] += 1
            showlegend = readId not in inLegend
            inLegend.add(readId)
            if thisIdentity >= lastIdentity:
                adjustedIdentity = lastIdentity - delta
            else:
                adjustedIdentity = thisIdentity
            lastIdentity = adjustedIdentity
            try:
                color = colors[readId]
            except KeyError:
                color = colors[readId] = next(colorIterator)
            scatters[readId].append(
                go.Scatter(
                    x=x,
                    y=[adjustedIdentity, adjustedIdentity],
                    name=readId,
                    text=f"{readId} {int(x[0])}-{int(x[1])} {thisIdentity:.4f}",
                    hoverinfo="text",
                    legendgroup=readId,
                    showlegend=showlegend,
                    mode=mode,
                    marker=marker,
                    line_color=color,
                )
            )

    # Change legend names to include the count of best windows the sequence is in.
    for readId, theseScatters in scatters.items():
        for scatter in theseScatters:
            scatter["name"] = f"{readId} ({bestCounts[readId]})"

    data = []

    # Add the scatters for each read in the order found in identity.
    for read in identity:
        data.extend(scatters[read.id])

    return data


def getAllLines(
    centers: list[float], identity: dict[DNARead, list[float]], args: argparse.Namespace
) -> list[go.Scatter]:
    """
    Plot identity lines for all sequences.

    @param centers: The x-axis centers of the genome windows.
    @param identity: A mapping from sequences (that have been matched against the
        reference) to a list of float identities (one for each window).
    @param title: The title of the plot.
    @param args: The command-line arguments.
    """
    data = []

    mode, marker = (
        ("lines+markers", dict(size=args.markerSize))
        if args.markerSize  # i.e., not 0 and not None.
        else ("lines", None)
    )

    for read, identities in identity.items():
        data.append(
            go.Scatter(
                x=centers,
                y=identities,
                mode=mode,
                marker=marker,
                name=read.id,
            )
        )

    return data


def plotPlotly(
    centers: list[float],
    identity: dict[DNARead, list[float]],
    title: str,
    args: argparse.Namespace,
) -> None:
    """
    Make a plotly plot.

    @param centers: The x-axis centers of the genome windows.
    @param identity: A mapping from sequences (that have been matched against the
        reference) to a list of float identities (one for each window).
    @param title: The title of the plot.
    @param args: The command-line arguments.
    """
    axisFontSize: int = 16
    titleFontSize: int = 18

    getData = getBestOnly if args.bestOnly else getAllLines
    data = getData(centers, identity, args)

    xaxis = {
        "title": "Genome location",
        "titlefont": {
            "size": axisFontSize,
        },
    }

    yaxis = {
        "title": "Nucleotide identity fraction",
        "titlefont": {
            "size": axisFontSize,
        },
    }

    layout = go.Layout(
        title={"text": title, "x": 0.5, "xanchor": "center"},
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            "size": titleFontSize,
        },
        hovermode="closest",
    )

    fig = go.Figure(data=data, layout=layout)

    if args.html:
        write_html(fig, args.html)

    if args.out:
        write_image(fig, args.out)


def main() -> None:
    """
    Read the FASTA, compute the windowed identity and then make either a matplotlib
    or plotly plot.
    """
    parser = makeParser()
    args = parser.parse_args()

    if not (args.out or args.html):
        sys.exit(
            "You must specify --out or --html to produce a static image or a "
            "plotly HTML image (or both).",
        )

    reads = tuple(
        DNARead(read.id, read.sequence) for read in parseFASTACommandLineOptions(args)
    )

    idsSeen: set[str] = set()
    for read in reads:
        if read.id in idsSeen:
            # We could be much more helpful here. E.g., compare the sequences, give the
            # index of the read, collect all errors and not fail on just the first.
            sys.exit(
                f"Read with id {read.id} seen more than once in the input.",
            )

    wi = WindowedIdentity(reads)

    reference, starts, identity = wi.getIdentity(
        referenceRegex=args.reference,
        window=args.window,
        minWindow=args.minWindow,
        jump=args.jump,
        includeRegex=args.include,
        sort=args.sort,
        strict=args.strict,
        startOffset=args.startOffset,
        stopOffset=args.stopOffset,
    )

    if args.dryRun:
        nToCompare = len(identity)
        s = "" if nToCompare == 1 else "s"
        print(f"The reference sequence is {reference.id!r}, ", end="")
        print(
            f"to which {nToCompare} sequence{s} would be compared:\n  "
            + "\n  ".join(sorted(read.id for read in identity))
        )
        sys.exit(0)

    # We will plot the identity fractions at the midpoints of the windows,
    # so we need the x-axis coordinates of the window centers.
    half = args.window / 2.0
    centers = [start + half for start in starts]

    title = args.title or (
        f"{args.window} nt (jump {args.jump}){' best' if args.bestOnly else ''} "
        f"windowed identity against {reference.id}"
    )

    if args.out and args.matplotlib:
        # Produce a static image using matplotlib.
        plotMatplotlib(centers, identity, title, args)

    if args.html or (args.out and not args.matplotlib):
        # Produce HTML (and maybe a static image) using plotly.
        plotPlotly(centers, identity, title, args)


if __name__ == "__main__":
    main()

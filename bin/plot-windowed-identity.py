#!/usr/bin/env python

import sys
import argparse
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.io import write_image, write_html

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
        required=True,
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
        "--matplotlib",
        action="store_true",
        help="Use matplotlib to make a static image. Otherwise, plotly will be used.",
    )

    parser.add_argument(
        "--y01",
        action="store_true",
        help="Set the y-axis limits to (0.0, 1.0)",
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

    return parser


def plotMatplotlib(
    reference: DNARead,
    centers: list[float],
    identity: dict[DNARead, list[float]],
    title: str,
    args: argparse.Namespace,
) -> None:
    """
    Make a matplotlib image.
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


def plotPlotly(
    reference: DNARead,
    centers: list[float],
    identity: dict[DNARead, list[float]],
    title: str,
    args: argparse.Namespace,
) -> None:
    """
    Make a plotly plot.
    """
    data = []
    axisFontSize: int = 16
    titleFontSize: int = 18

    for read, identities in identity.items():
        data.append(
            go.Scatter(
                x=centers,
                y=identities,
                # text=read.id,
                # hoverinfo="text",
                mode="lines",
                name=read.id,
            )
        )

    xaxis = {
        "title": "Genome location",
        "range": (0, len(reference)),
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
        f"{args.window} nt (jump {args.jump}) windowed identity against {reference.id}"
    )

    if args.out and args.matplotlib:
        # Produce a static image using matplotlib.
        plotMatplotlib(reference, centers, identity, title, args)

    if args.html or (args.out and not args.matplotlib):
        # Produce HTML (and maybe a static image) using plotly.
        plotPlotly(reference, centers, identity, title, args)


if __name__ == "__main__":
    main()

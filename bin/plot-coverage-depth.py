#!/usr/bin/env python

import argparse
import math
import sys
from itertools import cycle

import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from plotly.io import write_html, write_image

from dark.sam import UnknownReference, UnspecifiedReference, coverageDepth
from dark.utils import assignRegionRows, parseStartEndLabelColor

REGION_COLORS = px.colors.qualitative.Plotly


def niceTickStep(maxVal: float) -> int:
    """
    Return a round tick-step size that gives roughly 5–8 ticks over [0, maxVal].
    """
    if maxVal <= 0:
        return 1
    exp = math.floor(math.log10(maxVal))
    step = 10**exp
    if maxVal / step > 8:
        step *= 2
    elif maxVal / step < 4:
        step //= 2
    return max(1, int(step))


def makeParser() -> argparse.ArgumentParser:
    """
    Make an argument parser.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot per-base coverage depth from a SAM/BAM file.",
    )

    parser.add_argument(
        "sam",
        metavar="SAM-OR-BAM-FILE",
        help="The SAM or BAM file to read.",
    )

    parser.add_argument(
        "--out",
        required=True,
        help=(
            "Output filename. The file suffix determines the format: "
            ".html for interactive HTML, or .png/.svg/.pdf/.jpg for a static image."
        ),
    )

    parser.add_argument(
        "--reference",
        help=(
            "The reference sequence name to plot coverage for. If not given, "
            "the reference is auto-detected from the aligned reads (this fails "
            "if reads map to more than one reference or no mapped reads are found)."
        ),
    )

    parser.add_argument(
        "--title",
        help=(
            "The plot title. You can embed {median} {mean} {sd} {min} and {max} "
            "to incorporate coverage statistics into the title. The --places "
            "argument determines the number of decimal places that will be used in "
            "the mean and standard deviation. {reference} will embed the name of "
            "the reference."
        ),
    )

    parser.add_argument("--xTitle", default="Genome position", help="The x-axis title.")

    parser.add_argument("--yTitle", default="Coverage depth", help="The y-axis title.")

    parser.add_argument("--legendTitle", default="Regions", help="The legend title.")

    parser.add_argument(
        "--places",
        type=int,
        default=2,
        help=(
            "Number of decimal places for floating-point coverage statistics "
            "(median, mean, sd) that can be substituted into the plot title. See the "
            "--title option for how to do this."
        ),
    )

    parser.add_argument(
        "--titleFontSize",
        type=int,
        default=16,
        help=("The font size for the plot title."),
    )

    parser.add_argument(
        "--axisFontSize",
        type=int,
        default=14,
        help=("The font size for axis labels."),
    )

    parser.add_argument(
        "--legendSquareSize",
        type=int,
        default=8,
        help=("The length of a legend colour square (for --region arguments, if any)."),
    )

    parser.add_argument(
        "--region",
        action="append",
        dest="regions",
        metavar="START:END[:LABEL[:COLOR]]",
        help=(
            "Highlight a genome region with a coloured bar below y=0. "
            "START and END should be 1-based inclusive offsets. LABEL is optional "
            "text shown in the legend (default: region-1, region-2, …). "
            "COLOR is any CSS/plotly colour name or hex code; if omitted a "
            "colour is assigned automatically. Overlapping regions will be stacked "
            "vertically. This argument may be given multiple times."
        ),
    )

    return parser


def main() -> None:
    """
    Read coverage depth from a SAM/BAM file and produce a plotly plot.
    """
    parser = makeParser()
    args = parser.parse_args()

    # Parse region arguments.
    rawRegions = []
    if args.regions:
        for regionStr in args.regions:
            try:
                rawRegions.append(parseStartEndLabelColor(regionStr))
            except ValueError as e:
                sys.exit(f"Invalid --region argument: {e}")

    try:
        reference, depth = coverageDepth(args.sam, args.reference or None)
    except (UnknownReference, UnspecifiedReference) as e:
        sys.exit(str(e))

    refLength = len(depth)
    positions = list(range(1, refLength + 1))
    depthArray = np.array(depth)
    maxDepth = int(np.max(depthArray))

    places = args.places
    titleTemplate = args.title or f"Coverage depth: {reference}"
    title = titleTemplate.format(
        mean=f"{np.mean(depthArray):.{places}f}",
        median=f"{np.median(depthArray):.{places}f}",
        sd=f"{np.std(depthArray):.{places}f}",
        min=int(np.min(depthArray)),
        max=int(np.max(depthArray)),
        reference=reference,
    )

    # Fill in default labels and colors for regions.
    colorIter = cycle(REGION_COLORS)
    regions = []
    for i, (start, end, label, color) in enumerate(rawRegions):
        regions.append(
            (
                # Convert to 0-based.
                start - 1 or 0,
                end or refLength,
                label or f"region-{i + 1}",
                color or next(colorIter),
            )
        )

    # Assign overlapping regions to stacked rows.
    rowAssignments = assignRegionRows(regions) if regions else []
    nRows = max(rowAssignments) + 1 if rowAssignments else 0

    # Region bar geometry in data (y-axis) coordinates.
    # Row 0 sits just below y=0; successive rows stack downward.
    barHeight = max(1.0, maxDepth * 0.03)
    barGap = barHeight * 0.6

    def rowY(row: int) -> tuple[float, float]:
        """Return (y0, y1) for a region bar in the given row."""
        y1 = -(row * (barHeight + barGap) + barGap)
        y0 = y1 - barHeight
        return y0, y1

    yBottom = rowY(nRows - 1)[0] - barGap if nRows > 0 else 0.0
    yTop = maxDepth * 1.05

    # Y-axis ticks: only non-negative values.
    step = niceTickStep(maxDepth)
    tickVals = list(range(0, int(maxDepth) + step + 1, step))

    # Coverage area trace.
    coverageTrace = go.Scatter(
        x=positions,
        y=depth,
        mode="lines",
        fill="tozeroy",
        name="Coverage depth",
        showlegend=False,
        line=dict(color="steelblue", width=1),
        fillcolor="rgba(70, 130, 180, 0.3)",
    )

    # One phantom scatter per region for the legend, plus one shape per region.
    shapes = []
    legendTraces = []

    for (start, end, label, color), row in zip(regions, rowAssignments):
        y0, y1 = rowY(row)
        shapes.append(
            dict(
                type="rect",
                x0=start - 0.5,
                x1=end + 0.5,
                y0=y0,
                y1=y1,
                xref="x",
                yref="y",
                fillcolor=color,
                line_width=0,
                opacity=0.85,
            )
        )
        legendTraces.append(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=args.legendSquareSize, color=color, symbol="square"),
                # Add one to start to get a 1-based humanized start site.
                name=f"{label} ({start + 1}-{end})",
                showlegend=True,
            )
        )

    layout = go.Layout(
        title=dict(
            text=title,
            x=0.5,
            xanchor="center",
            font=dict(size=args.titleFontSize),
        ),
        xaxis=dict(
            title=dict(text=args.xTitle, font=dict(size=args.axisFontSize)),
            range=[0.5, refLength + 0.5],
        ),
        yaxis=dict(
            title=dict(text=args.yTitle, font=dict(size=args.axisFontSize)),
            range=[yBottom, yTop],
            tickmode="array",
            tickvals=tickVals,
            ticktext=list(map(str, tickVals)),
            zeroline=True,
            zerolinecolor="rgba(0,0,0,0.3)",
            zerolinewidth=1,
        ),
        shapes=shapes,
        hovermode="x unified",
        legend=dict(
            title=args.legendTitle,
            orientation="v",
        ),
    )

    fig = go.Figure(data=[coverageTrace] + legendTraces, layout=layout)

    if args.out.lower().endswith(".html"):
        write_html(fig, args.out)
    else:
        write_image(fig, args.out)


if __name__ == "__main__":
    main()

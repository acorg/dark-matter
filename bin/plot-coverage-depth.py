#!/usr/bin/env python

import argparse
import math
import sys
from itertools import cycle
from pathlib import Path

import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from plotly.io import write_html, write_image

from dark.sam import UnknownReference, UnspecifiedReference, coverageDepth

_REGION_COLORS = px.colors.qualitative.Plotly


def parseRegion(regionStr: str) -> tuple[int, int, str | None, str | None]:
    """
    Parse a --region argument of the form start:end[:label[:color]].

    @param regionStr: A C{str} region specification.
    @raise ValueError: If the region string is malformed.
    @return: A 4-tuple of (start, end, label, color) where start and end are
        1-based C{int} offsets (inclusive), label is a C{str} or C{None},
        and color is a C{str} or C{None}.
    """
    parts = regionStr.split(":", 3)
    if len(parts) < 2:
        raise ValueError(
            f"Region {regionStr!r} must contain at least a start and end "
            f"(e.g., '100:200')."
        )
    try:
        start = int(parts[0])
        end = int(parts[1])
    except ValueError:
        raise ValueError(
            f"Region {regionStr!r} start and end must be integers."
        )
    if start < 1:
        raise ValueError(
            f"Region {regionStr!r} start must be >= 1 (1-based coordinates)."
        )
    if end < start:
        raise ValueError(
            f"Region {regionStr!r} end ({end}) must be >= start ({start})."
        )
    label = parts[2] if len(parts) > 2 else None
    color = parts[3] if len(parts) > 3 else None
    return start, end, label, color


def assignRegionRows(regions: list) -> list[int]:
    """
    Assign each region to a row (0 = closest to y=0) using a greedy
    interval-scheduling algorithm so that no two regions in the same row
    overlap horizontally.

    @param regions: A list of (start, end, label, color) tuples.
    @return: A C{list} of 0-based C{int} row indices, one per region.
    """
    # Sort by start position so that non-overlapping regions are always
    # considered in left-to-right order, regardless of input order.
    indexed = sorted(enumerate(regions), key=lambda x: x[1][0])
    row_ends: list[int] = []
    assignments = [0] * len(regions)
    for orig_idx, (start, end, _, _) in indexed:
        placed = False
        for i, row_end in enumerate(row_ends):
            if start > row_end:
                row_ends[i] = end
                assignments[orig_idx] = i
                placed = True
                break
        if not placed:
            assignments[orig_idx] = len(row_ends)
            row_ends.append(end)
    return assignments


def niceTickStep(max_val: float) -> int:
    """
    Return a round tick-step size that gives roughly 5–8 ticks over [0, max_val].
    """
    if max_val <= 0:
        return 1
    exp = math.floor(math.log10(max_val))
    step = 10**exp
    if max_val / step > 8:
        step *= 2
    elif max_val / step < 4:
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
            "the reference is auto-detected from the aligned reads. This fails "
            "if reads map to more than one reference or if no mapped reads are found."
        ),
    )

    parser.add_argument(
        "--title",
        help="The plot title.",
    )

    parser.add_argument(
        "--xTitle",
        default="Genome position",
        help="The x-axis title.",
    )

    parser.add_argument(
        "--yTitle",
        default="Coverage depth",
        help="The y-axis title.",
    )

    parser.add_argument(
        "--legendTitle",
        default="Regions",
        help="The legend title.",
    )

    parser.add_argument(
        "--places",
        type=int,
        default=2,
        help=(
            "Number of decimal places for floating-point coverage statistics "
            "(mean, sd) substituted into the plot title."
        ),
    )

    parser.add_argument(
        "--region",
        action="append",
        dest="regions",
        metavar="START:END[:LABEL[:COLOR]]",
        help=(
            "Highlight a genomic region with a coloured bar below y=0. "
            "START and END are 1-based inclusive offsets. LABEL is optional "
            "text shown in the legend (default: region-1, region-2, …). "
            "COLOR is any CSS/plotly colour name or hex code; if omitted a "
            "colour is assigned automatically. Overlapping regions are stacked "
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
    raw_regions = []
    if args.regions:
        for regionStr in args.regions:
            try:
                raw_regions.append(parseRegion(regionStr))
            except ValueError as e:
                sys.exit(f"Invalid --region argument: {e}")

    # Compute coverage depth (auto-detects reference if not given).
    try:
        reference, depth = coverageDepth(args.sam, args.reference or None)
    except (UnknownReference, UnspecifiedReference) as e:
        sys.exit(str(e))

    refLength = len(depth)
    positions = list(range(1, refLength + 1))
    depth_arr = np.array(depth)
    max_depth = int(np.max(depth_arr))

    title_template = args.title or f"Coverage depth: {reference}"
    title = title_template.format(
        mean=f"{np.mean(depth_arr):.{args.places}f}",
        median=f"{np.median(depth_arr):.{args.places}f}",
        sd=f"{np.std(depth_arr):.{args.places}f}",
        min=int(np.min(depth_arr)),
        max=int(np.max(depth_arr)),
    )

    # Fill in default labels and colors for regions.
    color_iter = cycle(_REGION_COLORS)
    regions = []
    for i, (start, end, label, color) in enumerate(raw_regions):
        if label is None:
            label = f"region-{i + 1}"
        if color is None:
            color = next(color_iter)
        regions.append((start, end, label, color))

    # Assign overlapping regions to stacked rows.
    row_assignments = assignRegionRows(regions) if regions else []
    n_rows = max(row_assignments) + 1 if row_assignments else 0

    # Region bar geometry in data (y-axis) coordinates.
    # Row 0 sits just below y=0; successive rows stack downward.
    bar_height = max(1.0, max_depth * 0.03)
    bar_gap = bar_height * 0.6

    def row_y(row: int) -> tuple[float, float]:
        """Return (y0, y1) for a region bar in the given row."""
        y1 = -(row * (bar_height + bar_gap) + bar_gap)
        y0 = y1 - bar_height
        return y0, y1

    y_bottom = row_y(n_rows - 1)[0] - bar_gap if n_rows > 0 else 0.0
    y_top = max_depth * 1.05

    # Y-axis ticks: only non-negative values.
    step = niceTickStep(max_depth)
    tick_vals = list(range(0, int(max_depth) + step + 1, step))

    axisFontSize = 14
    titleFontSize = 16

    # Coverage area trace.
    coverage_trace = go.Scatter(
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
    legend_traces = []

    for (start, end, label, color), row in zip(regions, row_assignments):
        y0, y1 = row_y(row)
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
        legend_traces.append(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=12, color=color, symbol="square"),
                name=f"{label} ({start}-{end})",
                showlegend=True,
            )
        )

    layout = go.Layout(
        title=dict(
            text=title,
            x=0.5,
            xanchor="center",
            font=dict(size=titleFontSize),
        ),
        xaxis=dict(
            title=dict(text=args.xTitle, font=dict(size=axisFontSize)),
            range=[0.5, refLength + 0.5],
        ),
        yaxis=dict(
            title=dict(text=args.yTitle, font=dict(size=axisFontSize)),
            range=[y_bottom, y_top],
            tickmode="array",
            tickvals=tick_vals,
            ticktext=[str(v) for v in tick_vals],
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

    fig = go.Figure(data=[coverage_trace] + legend_traces, layout=layout)

    suffix = Path(args.out).suffix.lower()
    if suffix == ".html":
        write_html(fig, args.out)
    else:
        write_image(fig, args.out)


if __name__ == "__main__":
    main()

# Import and make available things that are convenient to have around in
# iPythonNotebook following 'from dark.ipynb import *'.

from .blast.alignments import BlastReadsAlignments
from .fasta import FastaReads
from .graphics import alignmentGraph, alignmentPanel, scoreGraph
from .html import (
    summarizeTitlesByCount,
    summarizeTitlesByLength,
    summarizeTitlesByMaxScore,
    summarizeTitlesByMedianScore,
    summarizeTitlesByTitle,
)
from .titles import TitlesAlignments, titleCounts

# Keep pyflakes quiet by pretending to make use of all our imports.
_ = (
    FastaReads,
    BlastReadsAlignments,
    titleCounts,
    TitlesAlignments,
    summarizeTitlesByLength,
    summarizeTitlesByMaxScore,
    summarizeTitlesByMedianScore,
    summarizeTitlesByCount,
    summarizeTitlesByTitle,
    alignmentGraph,
    alignmentPanel,
    scoreGraph,
)

__all__ = [
    "FastaReads",
    "BlastReadsAlignments",
    "titleCounts",
    "TitlesAlignments",
    "summarizeTitlesByLength",
    "summarizeTitlesByMaxScore",
    "summarizeTitlesByMedianScore",
    "summarizeTitlesByCount",
    "summarizeTitlesByTitle",
    "alignmentGraph",
    "alignmentPanel",
    "scoreGraph",
]

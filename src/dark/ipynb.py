# Import and make available things that are convenient to have around in
# iPythonNotebook following 'from dark.ipynb import *'.

from .fasta import FastaReads
from .blast.alignments import BlastReadsAlignments
from .titles import titleCounts, TitlesAlignments
from .html import (
    summarizeTitlesByLength,
    summarizeTitlesByMaxScore,
    summarizeTitlesByMedianScore,
    summarizeTitlesByCount,
    summarizeTitlesByTitle,
)
from .graphics import alignmentGraph, alignmentPanel, scoreGraph

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

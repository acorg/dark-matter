#!/usr/bin/env python

import sys
import os
import time
import argparse
import numpy as np
from sklearn.manifold import MDS
import plotly.express as px
import pandas as pd
from os.path import exists

from dark.sam import DistanceMatrix
from dark.utils import readLabels


# The following function is Hepatitis B virus (HBV) specific. So it really
# has no place being in the dark-matter repo. But it's convenient to have
# it here and I don't for the moment have a better (and certainly not
# simpler) place to put it and make it easily accessible (without adding
# another Python dependency), so here it sits, at least for now.
def HBVGenotypeKey(genotype):
    """
    Make a key for sorting HBV genotypes (including non-human genotypes that
    are the species name) for the plot legend.

    @param genotype: A C{str} genotype.
    @return: A 2-C{tuple} with an C{int} index and the passed C{str} genotype.
        This will ensure that sorting will be first on the index, then on the
        genotype name.
    """
    # The index is 0 for single-character genotypes (apart from ?, which
    # comes last), followed by 1 for genotypes with longer names
    # (Chimpanzee, Gibbon, etc), followed by 2 for unknown (?) genotypes.
    return (2 if genotype == "?" else (0 if len(genotype) == 1 else 1), genotype)


def getFigure(referenceIds, categories, dm, transform, args):
    """
    Make a 3D figure.

    @param referenceIds: A C{list} of C{str} reference ids.
    @param categories: A C{list} of C{str} reference id categories.
    @param dm: A C{DistanceMatrix} instance.
    @param transform: An C{MDS} transform.
    @param args: A C{Namespace} object returned by argparse, containing
        argument values passed on the command line.
    @return: A plotly express figure instance.
    """
    readCounts = [len(dm.scores[id_]) for id_ in referenceIds]

    if not (len(transform) == len(readCounts) == len(categories) == len(referenceIds)):
        raise ValueError(
            f"Unequal lengths: transform={len(transform)}, "
            f"readCounts={len(readCounts)}, categories={len(categories)}, "
            f"referenceIds={len(referenceIds)}."
        )

    if args.categorySortKey == "HBV genotype":
        categorySortKey = HBVGenotypeKey
    else:

        def categorySortKey(c):
            return c

    df = pd.DataFrame(
        {
            "x": transform[:, 0],
            "y": transform[:, 1],
            "Read count": readCounts,
            args.categoryName: categories,
            "Accession": referenceIds,
        }
    )

    categoryOrders = {args.categoryName: sorted(categories, key=categorySortKey)}

    if args.twoD:
        fig = px.scatter(
            df,
            x="x",
            y="y",
            color=args.categoryName,
            opacity=0.70,
            title="Reference map",
            hover_data=("Accession", "Read count"),
            category_orders=categoryOrders,
        )
    else:
        fig = px.scatter_3d(
            df,
            x="x",
            y="y",
            z="Read count",
            color=args.categoryName,
            opacity=0.70,
            title="Reference map",
            hover_data=("Accession",),
            category_orders=categoryOrders,
        )

    fig.update_traces(marker={"size": 4.2})

    fig.update_layout(
        margin=dict(l=0, r=0, b=0, t=0),
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )

    if args.logZ and not args.twoD:
        fig.update_layout(scene_zaxis_type="log")

    return fig


def getTransform(referenceIds, dm, args):
    """
    Get the MDS-transformed locations of the references.

    @param referenceIds: A C{list} of C{str} reference ids.
    @param dm: A C{DistanceMatrix} instance.
    @param args: A C{Namespace} object returned by argparse, containing
        argument values passed on the command line.
    @return: A numpy n x 2 array of x, y locations (where n is the number of
        reference ids).
    """
    if args.mdsFile and exists(args.mdsFile):
        if args.verbose:
            modificationTime = time.strftime(
                "%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(args.mdsFile))
            )
            print(
                f"Loading cached MDS transform file dated "
                f"{modificationTime} from {args.mdsFile!r}.",
                file=sys.stderr,
            )
        transform = np.load(args.mdsFile)
    else:
        distance = dm.matrix(
            referenceIds=referenceIds, metric=args.metric, similarity=False
        )
        mds = MDS(dissimilarity="precomputed", n_jobs=8)
        if args.verbose:
            print("Starting MDS fit.", file=sys.stderr)
        transform = mds.fit_transform(distance)
        if args.verbose:
            print("Finished MDS fit.", file=sys.stderr)

        # Cache the result if an MDS file name was given.
        if args.mdsFile:
            with open(args.mdsFile, "wb") as fp:
                np.save(fp, transform)

    return transform


def getDistanceMatrix(args):
    """
    Get the distance matrix for the references.

    @param args: A C{Namespace} object returned by argparse, containing
        argument values passed on the command line.
    @return: A C{DistanceMatrix} instance.
    """
    dm = DistanceMatrix()

    if args.samFile:
        for filename in args.samFile:
            dm.addFile(filename, scoreTag="AS")

        # Cache the distance matrix scores if a file name was given.
        if args.scoreFile:
            with open(args.scoreFile, "w") as fp:
                dm.save(fp)
    else:
        # The argparse group should ensure that --scoreFile was used if
        # --samFile was not.
        assert args.scoreFile

        if exists(args.scoreFile):
            with open(args.scoreFile) as fp:
                dm.load(fp)
        else:
            print("Score file {args.scoreFile!r} does not exist.", file=sys.stderr)
            sys.exit(1)

    return dm


def main(args):
    """
    Plot references according to their distances from each other, based on
    common read matches (and non-matches).

    @param args: A C{Namespace} object returned by argparse, containing
        argument values passed on the command line.
    """
    if not (args.samFile or args.scoreFile):
        print(
            "At least a SAM file or a score file name must be given.", file=sys.stderr
        )
        sys.exit(1)

    dm = getDistanceMatrix(args)

    if args.minMatchingReads is None:
        referenceIds = sorted(dm.scores)
    else:
        referenceIds = sorted(
            referenceId
            for (referenceId, reads) in dm.scores.items()
            if len(reads) >= args.minMatchingReads
        )
        if args.verbose:
            print(
                f"Found {len(referenceIds)} references with at least "
                f"{args.minMatchingReads} matching reads.",
                file=sys.stderr,
            )

    if args.categoryFile:
        with open(args.categoryFile) as fp:
            categoryDict = readLabels(fp)
        categories = [categoryDict.get(id_, id_) for id_ in referenceIds]
    else:
        categories = referenceIds

    # This plots a dendrogram - not needed / wanted at the moment.
    #
    # from matplotlib import pyplot as plt
    # from scipy.cluster.hierarchy import dendrogram, linkage
    # from scipy.spatial.distance import squareform
    #
    # condensed = squareform(distance)
    # linkage_ = linkage(condensed, method='average')
    # dendrogram(linkage_, orientation='left', labels=categories,
    #            leaf_font_size=2)
    # plt.show()

    transform = getTransform(referenceIds, dm, args)

    fig = getFigure(referenceIds, categories, dm, transform, args)

    fig.write_html(args.htmlFile)
    if args.verbose:
        print(f"Wrote {args.htmlFile!r}.", file=sys.stderr)


def makeParser():
    """
    Make a command-line argument parser.

    @return: An C{argparse.ArgumentParser} instance.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Plot references according to their distances from one "
            "another, based on common read matches (and "
            "non-matches)."
        )
    )

    parser.add_argument(
        "--samFile", action="append", help="The SAM file(s) to load. May be repeated."
    )

    parser.add_argument(
        "--scoreFile",
        help=(
            "A (JSON) file to read (JSON) the reference/read score matrix "
            "from. This should have been produced by an earlier call with a "
            "--scoreFile argument or else by saving the output of "
            "reference-read-scores-to-JSON.py  If the file does not exist, "
            "it will be created."
        ),
    )

    parser.add_argument(
        "--htmlFile", required=True, help="The HTML plotly file to write."
    )

    parser.add_argument(
        "--mdsFile", help="A JSON to cache the result of the MDS optimization."
    )

    parser.add_argument(
        "--categoryName",
        default="Genotype",
        help=(
            "The name of the categories to which reference sequences are "
            'assigned (e.g., "Genotype"). This name will appear at the '
            "top of the legend."
        ),
    )

    parser.add_argument(
        "--categorySortKey",
        choices=("HBV genotype",),
        help=(
            "The name of the function to use to sort categories (for the "
            "legend). If not given, categories will be sorted by name."
        ),
    )

    parser.add_argument(
        "--categoryFile",
        help=(
            "A file containing labels for reference sequences. The format "
            "is lines containing a reference name, a TAB, and a new name "
            "(to be shown in the plot). If no label file is given, or if "
            "there is no new name for a reference, the original reference "
            "names will be used."
        ),
    )

    parser.add_argument(
        "--minMatchingReads",
        type=int,
        help=(
            "The minimum number of reads that must match a reference for it "
            "to be included."
        ),
    )

    parser.add_argument(
        "--verbose", action="store_true", help="Print extra information."
    )

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--2d",
        "--2D",
        action="store_true",
        dest="twoD",
        help=(
            "Make the figure in two dimensions, instead of three. The third "
            "dimensions (reference matching read count) is omitted."
        ),
    )

    group.add_argument(
        "--logZ",
        action="store_true",
        help="Log the z (read count) axis (only valid if --2d is not used).",
    )

    parser.add_argument(
        "--metric",
        choices=("jaccard", "soergel"),
        default="soergel",
        help=("The distance metric to use."),
    )

    return parser


if __name__ == "__main__":
    main(makeParser().parse_args())

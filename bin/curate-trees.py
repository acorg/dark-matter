#!/usr/bin/env python

import sys
import argparse
from pathlib import Path
import dendropy
from warnings import warn


def parseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Curate phylogenetic trees (i.e., their specification files).",
    )

    parser.add_argument(
        "--treeFiles",
        nargs="+",
        help="The tree files to curate.",
    )

    parser.add_argument(
        "--outSuffix",
        default=".treefile-curated",
        help="The suffix to use on curated output files.",
    )

    parser.add_argument(
        "--dryRun",
        action="store_true",
        help="Just print what would be done.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite pre-existing result files.",
    )

    parser.add_argument(
        "--rootNode",
        help="Re-root the tree on the incoming branch to the tip node with this name.",
    )

    parser.add_argument(
        "--collapse",
        type=float,
        help=(
            "A floating-point support value. Branches with support less than this "
            "value will be collapsed (to form polytomies). The value to use will depend "
            "on what the program used to create the tree used to indicate the support. "
            "This is frequently the number of times the branch was supported in the "
            "bootstrap analysis. In this case, if you did 100 bootstraps, you would see "
            "integer support values from 0 to 100. If you did 1000 bootstraps, the "
            "values could go up to 1000. Other programs might write the support value "
            "as a fraction (of bootstraps in which the branch is supported), therefore "
            "ranging from 0.0 to 1.0. In other words, you need to know what's in your "
            "tree file to know for sure what value to use for this option. "
            "More difficult is to know what cut-off value actually makes sense in your "
            "analysis. If you're using iqtree2 with its ultrafast bootstrap (UFBoot), "
            "you should read the related FAQ item at "
            "http://www.iqtree.org/doc/Frequently-Asked-Questions"
        ),
    )

    parser.add_argument(
        "--discardLengths",
        action="store_true",
        help=(
            "When collapsing edges due to low support, the length of the edge being "
            "removed will be added to the length of the edges leading to the nodes of "
            "the edge-to-be-removed's child edges. This may be hard to follow, in "
            "which case you should draw yourself a little picture of part of a rooted "
            "tree and think about what happens to the nodes that become detached when "
            "you remove an edge (and the node that it leads to in the direction away "
            "from the root). If the length of the edge you are collapsing is "
            "important, because it came from a tree-making algorithm that sets branch "
            "lengths according to sequence change, the 'right' thing to do is to add "
            "the length of the to-be-deleted edge to the lengths of its descendant "
            "edges when you attach those nodes as children of the node at the root-end "
            "of the edge you are deleting. If for some reason you also want to discard "
            "the length of the edge you are collapsing, you can use this "
            "--discardLengths option. Make sure you really understand the implications "
            "of doing so, however, because discarding that length will cause the "
            "patristic distance of these child nodes to be shorter than it was "
            "estimated to be by the algorithm that made the tree and set the original "
            "branch lengths. And why would you want to second-guess the algorithm and "
            "mess with its conclusions?"
        ),
    )

    parser.add_argument(
        "--ladderize",
        choices=("ascending", "descending"),
        help="How to (optionally) ladderize the tree after re-rooting.",
    )

    parser.add_argument(
        "--format",
        default="newick",
        choices=("fasta", "phylip", "newick", "nexml", "nexus"),
        help="The output format.",
    )

    parser.add_argument(
        "--scale",
        default=1.0,
        type=float,
        help=(
            "Factor to scale the length of the branches leading to the descendants "
            "of the new root (if --rootNode is used). This can be used to reduce "
            "the length of very long branches that result from rerooting on a distant "
            "outgroup. If you use this option for an image in a publication, you of "
            "course must point out that you have done so!"
        ),
    )

    return parser.parse_args()


def main():
    args = parseArgs()

    if not (args.rootNode or args.ladderize or args.collapse):
        exit(
            "Nothing to do! Use some combination of --rootNode (for re-rooting), "
            "--ladderize, or --collapse if you want something done. Exiting."
        )

    for treefile in map(Path, args.treeFiles):
        curated = treefile.parent / f"{treefile.stem}{args.outSuffix}"
        if curated.exists() and not args.force:
            print(
                f"Output {str(curated)!r} exists. Skipping. Use --force to overwrite.",
                file=sys.stderr,
            )
            continue

        with open(treefile) as fp:
            tree = dendropy.Tree.get(
                file=fp, schema="newick", preserve_underscores=True
            )

        # Nodes must either have a taxon (for tips) or a 'label' (interpreted as the
        # support value on the incoming branch of internal nodes).
        for node in tree.preorder_node_iter():
            if node.label and node.taxon:
                # It's sufficiently weird to find a label on a tip branch that we should
                # wonder a) whether the code that created the support values was sane,
                # and b) whether these labels are in fact support values. If the value
                # is 100.0 there's no harm in ignoring it, but we issue a warning. If
                # we've been given a tree with branch labels that are numeric but not
                # support values, a taxon node with a label on it is a pretty good sign
                # that something is amiss.
                msg = (
                    f"Node {node!r} has a label and a taxon. It should have one or "
                    "the other, not both. The label is supposed to be a support "
                    "value, and tips (taxa) normally won't have one because they "
                    "always have 100% support."
                )
                if float(node.label) == 100.0:
                    warn(msg)
                else:
                    exit(msg)

        if args.collapse is not None:
            for node in tree.postorder_node_iter():
                if (
                    node.is_internal()
                    and node.label is not None
                    and float(node.label) < args.collapse
                ):
                    # Collapse this node's parent edge (i.e., the edge coming from the
                    # direction of the tree root) so its children become children of
                    # this node's parent.
                    node.edge.collapse(
                        adjust_collapsed_head_children_edge_lengths=not args.discardLengths
                    )

        # Re-root.
        if args.rootNode:
            node = tree.find_node_with_taxon_label(args.rootNode)
            assert node is not None
            tree.reroot_at_edge(
                node.edge,
                update_bipartitions=True,
                length1=node.edge.length / args.scale,
                length2=node.edge.length / args.scale,
            )

        if args.ladderize:
            # The reorder function is from dendropy v5.0.2 (released Sept 2024).
            tree.reorder()
            tree.ladderize(args.ladderize == "ascending")

        if args.dryRun:
            print(f"Would write curated tree file {str(curated)!r}.")
        else:
            tree.write(path=curated, schema=args.format)
            print(f"Wrote curated tree file {str(curated)!r}.")


if __name__ == "__main__":
    main()

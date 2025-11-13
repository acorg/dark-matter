#!/usr/bin/env python

import argparse
import sys
from typing import Optional

from ete4 import Tree

from dark.trees import ete4root


def parseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Curate a phylogenetic tree using the ete4 library and print the result."
        ),
    )

    parser.add_argument("treeFile", help="The tree file to curate.")

    parser.add_argument(
        "--parser",
        "--format",
        type=int,
        default=0,
        help=(
            "The Newick tree flavor. See "
            "https://etetoolkit.github.io/ete/reference/reference_parsers.html "
            "for information on the Newick variations ete4 understands. The output "
            "tree will be in the same format as the input."
        ),
    )

    parser.add_argument(
        "--rootNode",
        action="append",
        help=(
            "Re-root the tree on the branch entering the node that is the "
            "most-recent common ancestor of the tips with these names (may be "
            "repeated)."
        ),
    )

    parser.add_argument(
        "--scale",
        default=1.0,
        type=float,
        help=(
            "Factor to scale (i.e., divide) the length of the branches leading to the "
            "descendants of the new root (if --rootNode is used). This can be used to "
            "reduce the length of very long branches that result from rerooting on a "
            "distant outgroup. If you use this option for an image in a publication, "
            "you of course must point out that you have done so!"
        ),
    )

    parser.add_argument(
        "--detach",
        action="store_true",
        help=(
            "Detach the named root node from the tree after rooting on the edge coming "
            "into it. The newly added (by rooting) parent of the named node will be "
            "the (temporary) new root node. This option will first remove the named "
            "node as a child of that root, leaving it with just one child. The "
            "remaining child will then be made the new root.  This essentially like "
            "rooting on an outgroup taxon then completely removing the outgroup from "
            "the tree."
        ),
    )

    parser.add_argument(
        "--collapse",
        action="store_true",
        help="Collapse branches that are insufficiently supported.",
    )

    parser.add_argument(
        "--supportThreshold",
        type=int,
        help=(
            "The support level below which edges will be collapsed (only valid if "
            "--collapse is used)."
        ),
    )

    parser.add_argument(
        "--alrtThreshold",
        type=float,
        help=(
            "The alrt level below which edges will be collapsed (only valid if "
            "--collapse is used). Alrt values are added by iqtree2 if you use its "
            "-alrt option and can be used to check for well-supported branches."
        ),
    )

    parser.add_argument(
        "--removeSupport",
        action="store_true",
        help="Remove support values from internal nodes.",
    )

    parser.add_argument(
        "--removeAlrt",
        action="store_true",
        help="Remove alrt values from internal nodes.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print information about processing to standard error.",
    )

    return parser.parse_args()


def processSupport(tree: Tree, args: argparse.Namespace) -> None:
    """
    Potentially collapse and change the labels on internal nodes.
    """
    support: Optional[int]
    alrt: Optional[float]

    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            continue

        if node.name:
            try:
                alrts, supports = node.name.split("/")
            except ValueError:
                support = int(node.name)
                alrt = None
            else:
                alrt = float(alrts)
                support = int(supports)
        else:
            support = alrt = None

        if args.collapse:
            doCollapse = (
                args.supportThreshold is not None
                and support is not None
                and support < args.supportThreshold
            ) or (
                args.alrtThreshold is not None
                and alrt is not None
                and alrt < args.alrtThreshold
            )

            if doCollapse:
                if args.verbose:
                    print(f"Collapsing node with label {node.name}", file=sys.stderr)
                for child in node.children:
                    child.dist += node.dist

                node.delete()

        if args.removeSupport:
            if args.removeAlrt:
                name = ""
            else:
                name = "" if alrt is None else alrts
        else:
            if args.removeAlrt:
                name = supports
            else:
                name = node.name

        if node.name != name:
            if args.verbose:
                print(
                    f"Changing internal node label {node.name!r} -> {name!r}",
                    file=sys.stderr,
                )
            node.name = name


def main() -> None:
    args = parseArgs()
    tree = Tree(args.treeFile, parser=args.parser)

    processSupport(tree, args)

    if args.rootNode:
        tree = ete4root(tree, args.rootNode, args.detach, args.scale)

    print(tree.write(parser=args.parser))


if __name__ == "__main__":
    main()

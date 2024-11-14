#!/usr/bin/env python

import argparse
from ete3 import Tree


def parseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Curate a phylogenetic tree using the ete3 library and print the result."
        ),
    )

    parser.add_argument(
        "treeFile",
        help="The tree file to curate.",
    )

    parser.add_argument(
        "--format",
        type=int,
        default=0,
        help=(
            "The Newick tree flavor. See http://etetoolkit.org/docs/latest/tutorial/"
            "tutorial_trees.html#reading-and-writing-newick-trees for information on "
            "the Newick variations ete3 understands. The output tree will be in the "
            "same format as the input."
        ),
    )

    parser.add_argument(
        "--rootNode",
        required=True,
        help="Re-root the tree on the incoming branch to the tip node with this name.",
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

    return parser.parse_args()


def main():
    args = parseArgs()
    tree = Tree(args.treeFile, format=args.format)

    if roots := tree.search_nodes(name=args.rootNode):
        if len(roots) > 1:
            names = ", ".join(sorted(node.name for node in roots))
            exit(
                f"Found {len(roots)} tip nodes with names matching {args.rootNode!r}. "
                f"Found node names: {names}."
            )

        root = roots[0]
        tree.set_outgroup(root)

        if args.detach:
            parent = root.up
            assert len(parent.children) == 2
            root.detach()
            assert len(parent.children) == 1
            tree = parent.children[0]

        print(tree.write(format=args.format))
    else:
        exit(f"Could not find a tip named {args.rootNode!r}.")


if __name__ == "__main__":
    main()

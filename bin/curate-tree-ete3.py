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

    return parser.parse_args()


def main():
    args = parseArgs()
    tree = Tree(args.treeFile, format=args.format)

    if roots := tree.search_nodes(name=args.rootNode):
        if len(roots) > 1:
            exit(f"Could not find a tip named {args.rootNode!r}.")

        tree.set_outgroup(roots[0])
        print(tree.write(format=args.format))
    else:
        exit(f"Could not find a tip named {args.rootNode!r}.")


if __name__ == "__main__":
    main()

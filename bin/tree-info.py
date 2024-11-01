#!/usr/bin/env python

import argparse
import dendropy


def parseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Print information about a phylogenetic tree.",
    )

    parser.add_argument(
        "treeFile",
        help="The tree file to examine.",
    )

    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="The number of spaces to indent by when descending to children.",
    )

    parser.add_argument(
        "--lengths",
        action="store_true",
        help="Also show edge lengths.",
    )

    parser.add_argument(
        "--rootNode",
        help="Re-root the tree on the edge leading into the tip node with this name.",
    )

    return parser.parse_args()


def main() -> None:
    """
    Print the tree by calling the recursive node printer on the root node.
    """
    args = parseArgs()

    with open(args.treeFile) as fp:
        tree = dendropy.Tree.get(file=fp, schema="newick", preserve_underscores=True)

    if args.rootNode:
        node = tree.find_node_with_taxon_label(args.rootNode)
        assert node is not None
        tree.reroot_at_edge(node.edge, update_bipartitions=False)

    counts = dict.fromkeys(("tips", "internals"), 0)
    printNode(tree.seed_node, True, 0, args.indent, args.lengths, counts)

    totalNodes = counts["tips"] + counts["internals"]
    print(
        f"The tree has {totalNodes} nodes ({counts['internals']} "
        f"internal, {counts['tips']} tips).",
    )


def printNode(
    node: dendropy.Node,
    isRoot: bool,
    depth: int,
    indent: int,
    showLengths: bool,
    counts: dict[str, list[int, ...]],
) -> None:
    """
    Print a single node and (recursively) its descendants.

    @param node: The node to print.
    @param isRoot: True if this is the root of the tree.
    @param depth: The distance this node is from the root.
    @param indent: The number of spaces to indent the output by for each step away from
        the root.
    @param showLengths: If true, also print edge lengths.
    @param counts: A C{dict} with 'internals' and 'tips' keys used to accumulate the
        number of each node type during the recursive calls.
    """
    print((" " * depth * indent) + str(depth) + ": ", end="")

    length = f" {node.edge.length}" if showLengths else ""
    if node.is_internal():
        counts["internals"] += 1
        if node.label is None:
            if isRoot:
                print("Root")
            else:
                childLengths = [child.edge.length for child in node.child_nodes()]
                print(f"Int WARNING: no label{length}  Child lengths:", childLengths)
        else:
            print(f"Int {node.label}{length}")

        for child in node.child_nodes():
            printNode(child, False, depth + 1, indent, showLengths, counts)

    else:
        counts["tips"] += 1
        if node.taxon is None:
            print(f"Tip WARNING: no taxon{length}.")
        else:
            print(f"Tip {node.taxon.label}{length}")


if __name__ == "__main__":
    main()

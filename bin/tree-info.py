#!/usr/bin/env python

import argparse
import dendropy

# See https://jeetsukumaran.github.io/DendroPy/schemas/index.html for the available
# formats.
FORMATS = {"newick", "nexml", "nexus", "phylip"}


def parseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Print information phylogenetic tree nodes.",
    )

    parser.add_argument(
        "treeFile",
        help="The tree file to examine.",
    )

    parser.add_argument(
        "--format",
        default="newick",
        choices=sorted(FORMATS),
        help="The input tree format.",
    )

    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help=(
            "The number of spaces to increase node indentation by when descending to "
            "children."
        ),
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

    parser.add_argument(
        "--noTaxon",
        default="NO-TAXON!",
        help="What to print if a tip node has no taxon name.",
    )

    parser.add_argument(
        "--noLabel",
        default="NO-LABEL!",
        help="What to print if an internal node has no label.",
    )

    return parser.parse_args()


def main() -> None:
    """
    Print the tree by calling the recursive node printer on the root node.
    """
    args = parseArgs()

    kwargs = (
        dict(preserve_underscores=True) if args.format in {"newick", "nexus"} else {}
    )

    with open(args.treeFile) as fp:
        tree = dendropy.Tree.get(file=fp, schema=args.format, **kwargs)

    if args.rootNode:
        node = tree.find_node_with_taxon_label(args.rootNode)
        if not node:
            exit(f"Could not find tree node {args.rootNode!r} in {args.treeFile!r}.")
        tree.reroot_at_edge(node.edge, update_bipartitions=False)

    counts = dict.fromkeys(("tips", "internals"), 0)
    printNode(
        tree.seed_node, 0, args.indent, args.lengths, args.noTaxon, args.noLabel, counts
    )

    totalNodes = counts["tips"] + counts["internals"]
    print(
        f"The tree has {totalNodes} nodes ({counts['internals']} "
        f"internal, {counts['tips']} tips).",
    )


def printNode(
    node: dendropy.Node,
    depth: int,
    indent: int,
    showLengths: bool,
    noTaxon: str,
    noLabel: str,
    counts: dict[str, int],
) -> None:
    """
    Print a single node and (recursively) its descendants.

    @param node: The node to print.
    @param depth: The distance this node is from the root.
    @param indent: The number of spaces to indent the output by for each step away from
        the root.
    @param showLengths: If true, also print edge lengths.
    @param counts: A C{dict} with 'internals' and 'tips' keys used to accumulate the
        count of each node type during the recursive calls.
    """
    line = [f"{' ' * depth * indent} {depth}:"]

    length = node.edge.length if showLengths else ""
    if node.is_internal():
        counts["internals"] += 1
        children = node.child_nodes()
        line.append(f"Int {node.label or noLabel} {length}")
        if depth == 0:
            c = "child" if len(children) == 1 else 'children'
            line.append(f"(root has {len(children)} {c})")

        print(" ".join(line))

        for child in children:
            printNode(child, depth + 1, indent, showLengths, noTaxon, noLabel, counts)

    else:
        counts["tips"] += 1
        line.append(
            f"Tip {noTaxon if node.taxon is None else node.taxon.label} {length}"
        )
        print(" ".join(line))


if __name__ == "__main__":
    main()

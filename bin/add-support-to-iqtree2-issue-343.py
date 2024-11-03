#!/usr/bin/env python

import sys
import argparse
import dendropy

# See https://jeetsukumaran.github.io/DendroPy/schemas/index.html for the available
# formats.
FORMATS = {"newick", "nexml", "nexus", "phylip"}


def parseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Special-purpose script to add support labels to nodes in
trees produced by `iqtree2` when (if not run with `-keep-ident`) it adds nodes and tips
for identical sequences. iqtree2 currently does not put a support label on (the edge
leading to) the tip in the original processingor onto the nodes introduced by adding
tips for the identical sequences.  This is described at
https://github.com/iqtree/iqtree2/issues/343 The adjusted tree is written to standard
output in the same format as the input tree (as specified by --format argument).""",
    )

    parser.add_argument(
        "treeFile",
        help="The tree file to examine.",
    )

    parser.add_argument(
        "--support",
        default="100",
        help="The support label to add to unlabeled internal nodes",
    )

    parser.add_argument(
        "--format",
        default="newick",
        choices=sorted(FORMATS),
        help="The input (and resulting output) tree format.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print information about the number of nodes to which a label was added.",
    )

    return parser.parse_args()


def main() -> None:
    """
    Add a label (a support value) to internal nodes with no label and print the
    resultant tree to standard output.
    """
    args = parseArgs()

    kwargs = (
        dict(preserve_underscores=True) if args.format in {"newick", "nexus"} else {}
    )
    with open(args.treeFile) as fp:
        tree = dendropy.Tree.get(file=fp, schema=args.format, **kwargs)

    changes = 0

    for node in tree.preorder_internal_node_iter():
        if node.label is None:
            if node is not tree.seed_node:
                childLengths = [child.edge.length for child in node.child_nodes()]
                if childLengths != [0.0, 0.0]:
                    exit(
                        "Internal node found with edges to its children that are not "
                        f"all zero length (but instead, {', '.join(childLengths)}. "
                        "Are you sure this tree was produced by iqtree2?"
                    )
                node.label = args.support
                changes += 1

    if args.verbose:
        if changes:
            s = "" if changes == 1 else "s"
            print(
                f"Added label {args.support!r} to {changes} node{s}.",
                file=sys.stderr,
            )
        else:
            print(
                "No unlabeled nodes were found. The input and output trees will be "
                "identical.",
                file=sys.stderr,
            )

    tree.write(file=sys.stdout, schema=args.format)


if __name__ == "__main__":
    main()

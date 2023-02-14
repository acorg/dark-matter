#!/usr/bin/env python

import sys
from ete3 import Tree
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Print a Newick phylogenetic tree in ASCII.",
)

parser.add_argument(
    "--treeFile",
    metavar="FILENAME",
    default=sys.stdin,
    type=open,
    help="The Newick tree file.",
)

parser.add_argument(
    "--outgroupRegex",
    metavar="TAXON-REGEX",
    help=(
        "A regular expression (for taxa names). The lowest common parent of "
        "any matching taxa will be used as an outgroup."
    ),
)

parser.add_argument(
    "--format",
    metavar="N",
    type=int,
    help="The format of the Newick to be passed to ete3",
)

parser.add_argument(
    "--verbose",
    action="store_true",
    help=("Print information about the outgroup (if any) taxa to standard " "error"),
)

args = parser.parse_args()

tree = Tree(args.treeFile.read(), format=args.format)

if args.outgroupRegex:
    from re import compile

    regex = compile(args.outgroupRegex)
    taxa = [leaf.name for leaf in tree.iter_leaves() if regex.match(leaf.name)]

    if taxa:
        ca = tree.get_common_ancestor(taxa)
        if args.verbose:
            print("Taxa for outgroup:", taxa, file=sys.stderr)
            print("Common ancestor:", ca.name, file=sys.stderr)
            print("Common ancestor is tree:", tree == ca, file=sys.stderr)

        if len(taxa) == 1:
            tree.set_outgroup(tree & taxa[0])
        else:
            if ca == tree:
                tree.set_outgroup(tree.get_midpoint_outgroup())
            else:
                tree.set_outgroup(tree.get_common_ancestor(taxa))

print(tree.get_ascii())

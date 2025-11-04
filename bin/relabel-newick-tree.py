#!/usr/bin/env python

import argparse
import sys

from ete4 import Tree


def readNames(fp):
    """
    Read the existing / news names file.

    @param fp: A file pointer containing renaming information.

    @raise ValueError: If there is a duplicate old name given or if a line
        cannot be split into two TAB-separated fields.
    @return: A C{dict} of old to new names.
    """
    result = {}
    for line in fp:
        old, new = line.strip().split("\t")
        if old in result:
            raise ValueError(f"Sequence id {old!r} appears twice in the input.")
        else:
            result[old] = new

    return result


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Relabel a Newick tree according to a file of from/to names.",
)

parser.add_argument("newickFile", help="The Newick tree file.")

parser.add_argument(
    "--nameFile",
    help=(
        "The file of renamings. Each line should contain an existing (in "
        "the tree) sequence id, a TAB, and then a new sequence id. If not "
        "given, the renaming is read from standard input."
    ),
)

parser.add_argument(
    "--inputFormat",
    type=int,
    default=0,
    help=(
        "The parser argument to pass to the ete Tree class initializer to "
        "read the input Newick."
    ),
)

parser.add_argument(
    "--outputFormat",
    type=int,
    default=0,
    help="The parser argument to pass to the ete Tree writer for the resulting tree.",
)

parser.add_argument(
    "--quotedNames",
    action="store_true",
    help="If given, pass quoted_node_names=True to the ete Tree class initializer",
)

parser.add_argument(
    "--verbose",
    action="store_true",
    help="Print details of the renaming(s).",
)

args = parser.parse_args()

tree = Tree(
    args.newickFile, parser=args.inputFormat, quoted_node_names=args.quotedNames
)

if args.nameFile:
    with open(args.nameFile) as fp:
        names = readNames(fp)
else:
    names = readNames(sys.stdin)


# Relabel tree nodes that are mentioned in args.nameFile.
for node in tree.traverse():
    try:
        newName = names[node.name]
    except KeyError:
        pass
    else:
        if args.verbose:
            print(f"Renaming {node.name!r} to {newName!r}", file=sys.stderr)
        node.name = newName

# Print the tree, including the name of the root.
print("%s%s;" % (tree.write(parser=args.outputFormat)[:-1], tree.name))

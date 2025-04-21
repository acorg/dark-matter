import sys
import tempfile
import dendropy
from typing import TextIO


def removeAlrt(
    treeFile: str,
    supportThreshold: int,
    alrtThreshold: float,
    adjustSupportValues: bool,
    replacementSupportValue: int | str,
    expectedUnlabeledNodeCount: int,
    verbose: bool = False,
    verboseFp: TextIO = sys.stderr,
) -> str:
    """
    Read a Newick tree file produced by iqtree2 with both alrt and support
    values on nodes, drop the alrt values and optionally adjust the support
    values (based on the passed bootstrap support threshold and alrt
    threshold). Return the new tree as a (Newick) string.

    This function was written when using baltic
    (https://github.com/evogytis/baltic/tree/master) to draw trees.  Baltic
    cannot deal with node labels like 83.7/95 that are produced by iqtree2
    when you ask it for both bootstrap and alrt values (see
    http://www.iqtree.org/doc/Frequently-Asked-Questions). ete3 also cannot
    read those labels. But dendropy can. So we read the tree file with
    dendropy, examine the bootstrap support and alrt values, and return a
    tree with just (opionally adjusted) bootstrap support values.
    
    If C{adjustSupportValues} is True, the support on branches that have
    support values at least C{supportThreshold} but alrt values below
    C{alrtThreshold} will be set to the replacement support value.

    This allows you to put an artificially low support value onto such
    branches so that in display of the tree you can put a '*' onto branches
    with high support and nothing onto the ones with support below the
    threshold. I.e., this is a hack to take both support and alrt values into
    account in a context where you just want to display either a '*' or nothing
    at all on branches (the branches with the replacement support value will not
    have anything shown). Hopefully that makes sense!

    @param treeFile: The C{str} name of the Newick tree file to read.
    @param supportThreshold: The C{int} minimum support value needed to qualify for
        as high-support when C{adjustSupportValues} is true.
    @param alrtThreshold: The C{float} minimum alrt value needed to qualify for
        as high-support when C{adjustSupportValues} is true.
    @param adjustSupportValues: If C{True}, adjust the support on branches that have
        support values at least C{supportThreshold} but alrt values below
        C{alrtThreshold} to have a support value that is the passed support threshold
        minus one.
    @param replacementSupportValue: The C{int} or C{str} value to use as a label
        for branch whose support / alrt values are overall too low.
    @param expectedUnlabeledNodeCount: The number of unlabeled (internal) nodes
        expected. Pass -1 to allow any number of unlabeled internal nodes.
    @param verbose: If true, print information about nodes with good support but low
        alrt values to the C{verboseFp} output stream.
    @param verboseFp: The open output stream to write to if C{verbose} is true.

    """
    replacementSupportLabel = str(
        supportThreshold - 1
        if replacementSupportValue is None else replacementSupportValue
    )

    with open(treeFile) as fp:
        tree = dendropy.Tree.get(file=fp, schema="newick", preserve_underscores=True)

    unsupportedDueToLowAlrt = unlabeledNodeCount = nodeCount = 0

    for nodeCount, node in enumerate(
        tree.preorder_internal_node_iter(exclude_seed_node=False), start=1
    ):
        try:
            alrts, supports = node.label.split("/")
        except ValueError:
            raise ValueError(
                f"WARNING: node value ({node.label!r}) on node {node!r} could "
                "not be split on '/'.",
            )
        except AttributeError:
            # This can happen because the tree has been rooted and one edge
            # leading away from the inserted root node leads to a node with
            # no support label (this always happens when rooting on an edge
            # that leads to a tip). Or because a root node was created with
            # no label.
            assert node.label is None
            unlabeledNodeCount += 1
            continue
        else:
            node.label = supports
            if adjustSupportValues:
                support = int(supports)
                alrt = float(alrts)
                if support >= supportThreshold and alrt < alrtThreshold:
                    # When drawn, the edge coming to this node should not
                    # have a star, because its combined support / alrt values
                    # aren't high enough.  So we set its support label to the
                    # special value.
                    node.label = replacementSupportLabel
                    unsupportedDueToLowAlrt += 1

    if (
        expectedUnlabeledNodeCount != -1
        and unlabeledNodeCount != expectedUnlabeledNodeCount
    ):
        raise ValueError(
            f"Found {unlabeledNodeCount} internal nodes with no label "
            f"in {treeFile!r}, but {expectedUnlabeledNodeCount} were expected."
        )

    # Sanity check that there are now no node labels with a slash in them.
    for node in tree.preorder_internal_node_iter(exclude_seed_node=True):
        assert node.label is None or "/" not in node.label

    if verbose:
        print(
            f"Of {nodeCount} internal nodes, {unsupportedDueToLowAlrt} had "
            f"an alrt value that did not meet the {alrtThreshold} threshold, "
            "despite having a sufficiently high bootstrap support "
            f"({supportThreshold}) to otherwise qualify as well-supported.",
            file=verboseFp,
        )

    with tempfile.NamedTemporaryFile() as new:
        tree.write(path=new.name, schema="newick")
        result = open(new.name).read().strip()

    return result

#!/usr/bin/env python

import argparse

from dark.trees import removeAlrt


def makeParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            """Read a Newick tree file produced by iqtree2 with both alrt
and support values on nodes, drop the alrt values and optionally
adjust the support values (based on the passed bootstrap support
threshold and alrt threshold). Print the new tree as a (Newick) string."""
        ),
    )

    parser.add_argument(
        "treefile",
        help="The Newick tree file to read.",
    )

    parser.add_argument(
        "--supportThreshold",
        default=95,
        type=int,
        help="The integer minimum UF boot support value.",
    )

    parser.add_argument(
        "--alrtThreshold",
        default=80.0,
        type=float,
        help="The float minimum alrt value.",
    )

    parser.add_argument(
        "--adjustSupportValues",
        action="store_true",
        help=(
            "The support on branches that have support values at least "
            "--supportThreshold but alrt values below --alrtThreshold "
            "will be set to the --replacementSupportValue value."
        ),
    )

    parser.add_argument(
        "--replacementSupportValue",
        help=(
            "The string value to use as a label on branches whose support "
            "values are adjusted."
        ),
    )

    parser.add_argument(
        "--expectedUnlabeledNodeCount",
        type=int,
        help="The number of unlabeled nodes (in the input file) to allow.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "Print information about how many branches have good support "
            "but too-low alrt values."
        ),
    )

    return parser


def main() -> None:
    args = makeParser().parse_args()

    print(
        removeAlrt(
            args.treefile,
            args.supportThreshold,
            args.alrtThreshold,
            args.adjustSupportValues,
            args.replacementSupportValue,
            args.expectedUnlabeledNodeCount,
            args.verbose,
        )
    )


if __name__ == "__main__":
    main()

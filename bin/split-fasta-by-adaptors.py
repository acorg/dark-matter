#!/usr/bin/env python

from Bio import SeqIO
from collections import defaultdict
import sys
from dark.distance import levenshtein
from math import log10, ceil

# The name of the unknown adaptor.
UNKNOWN = "UNKNOWN"


def splitFASTAByAdaptor(
    knownAdaptors,
    adaptorLen,
    adaptorOffset,
    maximumDistance,
    outputPrefix,
    dryRun,
    verbose,
):
    """
    @param knownAdaptors: A C{set} of expected adaptor sequences.
    @param adaptorLen: The C{int} length of each adaptor sequence.
    @param adaptorOffset: The zero-based C{int} offset of the adaptor in
        each sequence.
    @param maximumDistance: The maximum distance an unknown adaptor will be
        mapped to in an attempt to find its nearest known adaptor.
    @param outputPrefix: A C{str} prefix that should be used in the file names
        that are written out.
    @param dryRun: A C{bool}, if C{True} only print what would be done, don't
        create any new FASTA files.
    @param verbose: A C{bool}, if C{True} output additional information about
        adaptor classes found and assigned.
    """
    adaptors = defaultdict(int)
    unknowns = 0
    classes = dict(zip(knownAdaptors, knownAdaptors))
    reads = []

    for count, seq in enumerate(SeqIO.parse(sys.stdin, "fasta"), start=1):
        reads.append(seq)
        adaptor = str(seq.seq)[adaptorOffset:][:adaptorLen].upper()
        adaptors[adaptor] += 1

    order = sorted(adaptors, key=lambda adaptor: adaptors[adaptor], reverse=True)

    for adaptor in order:
        if adaptor in knownAdaptors:
            if verbose:
                print("%s: %s. Known adaptor" % (adaptor, adaptors[adaptor]))
        else:
            distances = sorted(
                (levenshtein(adaptor, known), known) for known in knownAdaptors
            )
            # Treat the read as unclassifiable if it's too far from its
            # nearest neighbor or if its nearest neighbor is ambiguous.
            nearest = distances[0][0]
            if nearest > maximumDistance or (
                len(knownAdaptors) > 1 and nearest == distances[1][0]
            ):
                unknowns += 1
                classes[adaptor] = UNKNOWN
                if verbose:
                    print(
                        "%s: %s. Unknown, distances %r"
                        % (adaptor, adaptors[adaptor], [d[0] for d in distances])
                    )
            else:
                correctedAdaptor = distances[0][1]
                classes[adaptor] = correctedAdaptor
                if verbose:
                    print(
                        "%s: %s. Assigned to class %s, at dist %d"
                        % (
                            adaptor,
                            adaptors[adaptor],
                            correctedAdaptor,
                            distances[0][0],
                        )
                    )

    readGroups = defaultdict(list)

    # Collect reads into classes.
    for read in reads:
        adaptor = str(read.seq)[adaptorOffset:][:adaptorLen].upper()
        readGroups[classes[adaptor]].append(read[adaptorOffset + adaptorLen :])

    # Calculate the number of digits in the size of the biggest read group
    # so we can nicely align the output.
    width = int(ceil(log10(max(len(group) for group in readGroups.values()))))

    # The width of the count of files we'll write, so file names have zero
    # padded numeric prefixes.
    filesWidth = int(ceil(log10(len(readGroups))))

    # Write out the FASTA files for each adaptor class (this includes the
    # unclassifiable reads if any unknown adaptors were found).
    for count, adaptor in enumerate(sorted(readGroups), start=1):
        reads = readGroups[adaptor]
        filename = "%s%0*d-%s.fasta" % (outputPrefix, filesWidth, count, adaptor)
        description = (
            "unrecognized adaptors" if adaptor == UNKNOWN else "adaptor %s" % adaptor
        )
        if dryRun:
            print(
                "Would write %*d sequences for %s to %s"
                % (width, len(reads), description, filename)
            )
        else:
            with open(filename, "w") as fp:
                SeqIO.write(reads, fp, "fasta")
            print(
                "Wrote %*d sequences for %s to %s"
                % (width, len(reads), description, filename)
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "From a set of known adaptors and a FASTA on stdin, "
            "write FASTA files split by adaptor"
        )
    )

    parser.add_argument(
        "adaptors",
        nargs="+",
        metavar="adaptor",
        help="the set of adaptors that were used in sequencing",
    )

    parser.add_argument(
        "--adaptor-offset",
        metavar="N",
        type=int,
        default=0,
        dest="adaptorOffset",
        help="the offset of the adaptor in each sequence",
    )

    parser.add_argument(
        "--maximum-distance",
        metavar="dist",
        type=int,
        default=2,
        dest="maximumDistance",
        help="The maximum distance an unknown adaptor will be "
        "mapped to in an attempt to find its nearest known adaptor.",
    )

    parser.add_argument(
        "--outputPrefix",
        metavar="prefix",
        default="MID-",
        help="the prefix to use for newly created FASTA files",
    )

    parser.add_argument(
        "--dry-run",
        type=bool,
        default=False,
        dest="dryRun",
        help="If True, do not write new FASTA files, just show what would be " "done",
    )

    parser.add_argument(
        "--verbose",
        type=bool,
        default=False,
        help="If True, print information about adaptor classes",
    )

    args = parser.parse_args()
    adaptorLen = len(args.adaptors[0])
    adaptors = set(args.adaptors)

    # Check all adaptors are the same length.
    if any((len(adaptor) != adaptorLen) for adaptor in adaptors):
        print("All adaptors must be the same length.", file=sys.stderr)
        sys.exit(1)

    splitFASTAByAdaptor(
        adaptors,
        adaptorLen,
        args.adaptorOffset,
        args.maximumDistance,
        args.outputPrefix,
        args.dryRun,
        args.verbose,
    )

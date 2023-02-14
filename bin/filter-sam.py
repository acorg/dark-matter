#!/usr/bin/env python

import sys
from pysam import AlignmentFile

from dark.filter import (
    addFASTAFilteringCommandLineOptions,
    parseFASTAFilteringCommandLineOptions,
)
from dark.reads import Reads
from dark.sam import samfile, SAMFilter


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Given a SAM/BAM file and a set of filtering criteria "
            "write filtered SAM/BAM to stdout."
        ),
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="If True, do not print the final summary.",
    )

    parser.add_argument(
        "--bam",
        action="store_const",
        const="b",
        default="",
        help="If given, write (gzip compressed) BAM output.",
    )

    parser.add_argument(
        "--checkResultCount",
        type=int,
        help=(
            "The number of alignments expected in the output. If this "
            "number is not seen, the script exits with status 1 (and an "
            "error message is also printed, unless --quiet was used)."
        ),
    )

    addFASTAFilteringCommandLineOptions(parser)
    SAMFilter.addFilteringOptions(parser)

    args = parser.parse_args()
    reads = parseFASTAFilteringCommandLineOptions(args, Reads())
    samFilter = SAMFilter.parseFilteringOptions(
        args, reads.filterRead, storeQueryIds=True
    )

    # The following 'if' has a False in it to make it always fail. That's
    # because pysam issue 716 (see below) did not fix the problem as I had
    # hoped. Instead it throws an error if you pass a header that has a
    # modified SQ key with reference ids and there's a difference it
    # doesn't like. It's always safe to use the 'else' below, with the
    # slight downside being that its header will mention all sequence ids,
    # even if you only want a lesser number (via --referenceId). I'm
    # leaving the code here because this is how you would do it, and it
    # might be possible to just copy the 'header' dict below and further
    # adjust it to avoid the pysam error.
    if False and samFilter.referenceIds:
        # Make a header that only includes the wanted reference ids (if
        # any).
        print(
            "If you get a segmentation violation, this is (I think) "
            "due to the following pysam issue "
            "https://github.com/pysam-developers/pysam/issues/716 It is "
            "safer to run without using --referenceId.",
            file=sys.stderr,
        )
        with samfile(args.samfile) as sam:
            header = sam.header.as_dict()
        referenceIds = samFilter.referenceIds
        sequences = [
            sequence for sequence in header["SQ"] if sequence["SN"] in referenceIds
        ]
        referenceIdsFound = set(sequence["SN"] for sequence in sequences)
        notFound = referenceIds - referenceIdsFound
        if notFound:
            if len(notFound) == len(referenceIds):
                raise ValueError(
                    "None of the provided --referenceId ids were found in %s"
                    % args.samfile
                )
            else:
                # Just warn about the missing references.
                print(
                    "WARNING: reference%s (%s) given with --referenceId not "
                    "found in %s."
                    % (
                        "" if len(notFound) == 1 else "s",
                        ", ".join(sorted(notFound)),
                        args.samfile,
                    ),
                    file=sys.stderr,
                )

        header["SQ"] = sequences
        out = AlignmentFile(sys.stdout, mode="w" + args.bam, header=header)
    else:
        with samfile(args.samfile) as sam:
            out = AlignmentFile(sys.stdout, mode="w" + args.bam, template=sam)

    save = out.write
    kept = 0
    for kept, alignment in enumerate(samFilter.alignments(), start=1):
        save(alignment)

    out.close()

    if not args.quiet:
        total = samFilter.alignmentCount
        print(
            "Read %d alignment%s, kept %d (%.2f%%)."
            % (
                total,
                "" if total == 1 else "s",
                kept,
                0.0 if total == 0 else kept / total * 100.0,
            ),
            file=sys.stderr,
        )

    if args.checkResultCount is not None and kept != args.checkResultCount:
        if not args.quiet:
            print(
                "Did not write the expected %d alignment%s (wrote %d)."
                % (
                    args.checkResultCount,
                    "" if args.checkResultCount == 1 else "s",
                    kept,
                ),
                file=sys.stderr,
            )
        sys.exit(1)

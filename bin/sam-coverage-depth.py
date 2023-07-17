#!/usr/bin/env python

import sys
import argparse
from numpy import std
from collections import defaultdict
from json import dumps

from dark.filter import (
    addFASTAFilteringCommandLineOptions,
    parseFASTAFilteringCommandLineOptions,
)
from dark.reads import Reads
from dark.sam import samfile, SAMFilter, samReferences, UnknownReference
from dark.utils import baseCountsToStr, pct

from gb2seq import Gb2SeqError
from gb2seq.alignment import Gb2Alignment, addAlignerOption
from gb2seq.features import Features, addFeatureOptions


def makeParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Print SAM/BAM file coverage statistics by offset. "
            "Output lines show the offset."
        ),
    )

    addFASTAFilteringCommandLineOptions(parser)
    SAMFilter.addFilteringOptions(parser, samfileIsPositional=True)

    parser.add_argument(
        "--noOffsets",
        action="store_true",
        help="Do not print per-offset details of base counts.",
    )

    parser.add_argument(
        "--noStats",
        action="store_true",
        help="Do not print final average and standard deviation statistics.",
    )

    parser.add_argument(
        "--noFilter",
        action="store_true",
        help=(
            "Do not use our SAM filtering. Note that if you give this option, "
            "any filtering option (other than --referenceId) you also specify "
            "that is provided by the SAMFilter.addFilteringOptions will be "
            "silently ignored!"
        ),
    )

    addAlignerOption(parser)
    addFeatureOptions(
        parser,
        referenceHelpInfo=(
            " If you use this option, output will be printed as one "
            "JSON object per line, giving information not just about "
            "base coverage at each site but also about features in "
            "the reference and genome."
        ),
    )

    return parser


def printSiteInfo(
    alignment, site, includeUntranslated, includeGenome, bases, baseCount
):
    """
    Report what's found at a site for a given genome (or report insufficient
    coverage to standard error).

    @param alignment: A C{Gb2Alignment} instance.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @param includeGenome: If C{True}, include information about the genome
        (not just the reference).
    @param bases: A C{dict} mapping C{str} base names to C{int} counts.
    @param baseCount: The C{int} number of bases at the site (i.e., the
        coverage).
    """
    try:
        offsetInfo = alignment.offsetInfo(
            site,
            relativeToFeature=False,
            aa=False,
            featureName=None,
            includeUntranslated=includeUntranslated,
            allowAmbiguous=True,
        )
    except Gb2SeqError as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    if offsetInfo is None:
        print("No coverage")
        return

    if not includeGenome:
        del offsetInfo["genome"]

    offsetInfo["featureNames"] = sorted(offsetInfo["featureNames"])

    # Add an "nt" entry to the reference and genome and delete the id.
    gappedOffset = alignment.gappedOffsets[site]
    for what, gappedSequence in zip(
        ("reference", "genome"),
        (alignment.referenceAligned.sequence, alignment.genomeAligned.sequence),
    ):
        info = offsetInfo[what]
        try:
            info["nt"] = gappedSequence[gappedOffset]
        except IndexError:
            error = (
                f"Index error determining nt at site {site} of {what}: "
                f"{gappedOffset=}. Alignment has length {len(gappedSequence)}. "
                f"Ref has length {len(alignment.features.reference)}."
            )
            # An index error was occurring before we started passing start,
            # stop, and truncate=True to pysam's pileup method, which was
            # returning an offset bigger than the reference genome
            # length. See
            # https://github.com/VirologyCharite/bih-pipeline/issues/199
            # (assuming you have access to that repo). This should now be
            # fixed. Here, just print a message to standard error, put the
            # error into the JSON that is printed, and keep going.
            print(error, file=sys.stderr)
            print(offsetInfo, file=sys.stderr)
            info["IndexError"] = error
        del offsetInfo[what]["id"]

    offsetInfo["bases"] = bases
    offsetInfo["coverage"] = baseCount

    print(dumps(offsetInfo, indent=None, sort_keys=True))


def main():
    parser = makeParser()
    args = parser.parse_args()

    if args.noOffsets and args.noStats:
        print(
            "You have used both --noOffsets and --noStats, so there is no " "output!",
            file=sys.stderr,
        )
        sys.exit(1)

    # We don't have a file of reads, we just want a read filter that we can use
    # to filter the SAM file query sequences and to get reference lengths from.
    reads = parseFASTAFilteringCommandLineOptions(args, Reads())
    samFilter = SAMFilter.parseFilteringOptions(args, reads.filterRead)

    printOffsets = not args.noOffsets
    printStats = not args.noStats

    if samFilter.referenceIds and len(samFilter.referenceIds) > 1:
        print(
            "Only one reference id can be given. To calculate coverage for more "
            "than one reference, run this script multiple times.",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        referenceLengths = samFilter.referenceLengths()
    except UnknownReference:
        referenceId = samFilter.referenceIds.pop()
        referenceIds = samReferences(args.samfile)
        print(
            "Reference %r does not appear in SAM file %s. Known "
            "references are: %s."
            % (referenceId, args.samfile, ", ".join(sorted(referenceIds))),
            file=sys.stderr,
        )
        sys.exit(1)

    if args.noFilter:
        # Do not do our custom SAM filtering.
        def filterRead(read):
            return not (read.is_del or read.is_refskip)

    else:

        def filterRead(read):
            return not (read.is_del or read.is_refskip) and samFilter.filterAlignment(
                read.alignment
            )

    if printStats:
        counts = []

    if args.reference or args.sars2:
        features = Features(
            args.reference,
            sars2=args.sars2,
            addUnannotatedRegions=args.addUnannotatedRegions,
        )

        alignment = Gb2Alignment(features.reference, features, aligner=args.aligner)
    else:
        alignment = None

    with samfile(args.samfile) as sam:
        if samFilter.referenceIds:
            # No need to check if the given reference id is in referenceLengths
            # because the samFilter.referenceLengths call above catches that.
            referenceId = samFilter.referenceIds.pop()
        else:
            if len(referenceLengths) == 1:
                referenceId = list(referenceLengths)[0]
            else:
                print(
                    "SAM file %r contains %d references (%s). Only one "
                    "reference id can be analyzed at a time. Please use "
                    "--referenceId to specify the one you want examined."
                    % (
                        args.samfile,
                        len(referenceLengths),
                        ", ".join(sorted(referenceLengths)),
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)

        for column in sam.pileup(
            reference=referenceId,
            start=0,
            stop=referenceLengths[referenceId],
            truncate=True,
        ):
            bases = defaultdict(int)
            for read in column.pileups:
                if filterRead(read):
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base] += 1

            baseCount = sum(bases.values())

            if printStats:
                counts.append(baseCount)

            if printOffsets:
                if alignment:
                    printSiteInfo(
                        alignment,
                        column.reference_pos,
                        args.addUnannotatedRegions,
                        True,
                        bases,
                        baseCount,
                    )
                else:
                    print(
                        f"{column.reference_pos + 1} {baseCount} "
                        f"{baseCountsToStr(bases)}"
                    )

    if printStats:
        referenceLength = referenceLengths[referenceId]
        print("Reference id: %s" % referenceId)
        print("Reference length: %d" % referenceLength)
        print("Bases covered: %s" % pct(len(counts), referenceLength))
        print(
            "Min coverage depth: %d"
            % (0 if len(counts) < referenceLength else min(counts))
        )
        if counts:
            # Don't use Python3 default= option on max. Trying to keep Python 2
            # compatibility.
            print("Max coverage depth: %d" % max(counts))
        else:
            print("Max coverage depth: 0")
        print("Mean coverage depth: %.3f" % (sum(counts) / referenceLength))
        if counts:
            print("Coverage depth s.d.: %.3f" % std(counts))


if __name__ == "__main__":
    main()

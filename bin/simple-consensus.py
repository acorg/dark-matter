#!/usr/bin/env python

import sys
import argparse

from dark import __version__
from dark.consensus import consensusFromBAM, ConsensusError
from dark.sam import SamError


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Make a consensus sequence.",
    )

    parser.add_argument(
        "--bamFilename",
        required=True,
        help="The BAM file from which the consensus should be called.",
    )

    parser.add_argument(
        "--version", action="store_true", help="Print the version number and exit"
    )

    parser.add_argument(
        "--bamId",
        metavar="ID",
        help=(
            "The BAM file reference name indicating which aligned "
            "reads to make a consensus from. If not given, will be inferred "
            "from the BAM file header."
        ),
    )

    parser.add_argument(
        "--referenceFasta", metavar="FILE", help="The reference FASTA file."
    )

    parser.add_argument(
        "--fastaId",
        metavar="ID",
        help=(
            "The id of the sequence in --referenceFasta to use as a "
            "reference. Only considered if --referenceFasta is used. If not "
            "given and --referenceFasta is, the reference id will be "
            "inferred from reference names in the BAM header, or will be "
            "taken as the id of the first sequence in --referenceFasta."
        ),
    )

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--consensusId",
        metavar="ID",
        help=(
            "The id to use in the consensus sequence in the output FASTA. "
            "If not given, the reference sequence id will be used."
        ),
    )

    group.add_argument(
        "--idLambda",
        metavar='"lambda id: ..."',
        help=(
            "A one-argument function taking and returning a read id. "
            "This can be used to set the id of the consensus sequence based "
            "on the id of the reference sequence. The function will be "
            "called with the id of the BAM reference sequence. E.g., "
            "--idLambda \"lambda id: id.split('_')[0]\" or "
            "--idLambda \"lambda id: id[:10] + '-consensus'\"."
        ),
    )

    parser.add_argument(
        "--minCoverage",
        default=0,
        type=int,
        metavar="N",
        help=(
            "The minimum number of reads that must cover a site for a "
            "consensus base to be called. If zero reads cover a site, the "
            "--noCoverage value is used or if the number is greater than "
            "zero but less than minCoverage, the --lowCoverage value is "
            "used."
        ),
    )

    parser.add_argument(
        "--lowCoverage",
        default="n",
        metavar="N",
        help=(
            "What to do when some reads cover a site, but fewer than "
            'the --minCoverage value. Either set this to "reference" (to '
            "use the base from the reference sequence)or a single character "
            '(e.g., "N").'
        ),
    )

    parser.add_argument(
        "--noCoverage",
        default="N",
        help=(
            "What to do when no reads cover a site, Either set this to "
            '"reference" (to use the base from the reference sequence), or '
            'a single character (e.g., "N").'
        ),
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.6,
        metavar="N",
        help=(
            "The frequency threshold when calling the consensus. If the "
            "frequency of the most-common nucleotide at a site meets this "
            "threshold, that nucleotide will be called. Otherwise, an "
            "ambiguous nucleotide code will be produced, based on the "
            "smallest set of most-frequent nucleotides whose summed "
            "frequencies meet the threshold. If the frequency of the "
            "nucleotide that causes the threshold to be reached is the "
            "same as that of other nucleotides, all such nucleotides will "
            "be included in the ambiguous code."
        ),
    )

    parser.add_argument(
        "--ignoreQuality",
        action="store_true",
        help=(
            "Ignore FASTQ quality scores. All bases of all reads will be "
            "considered equally reliable."
        ),
    )

    parser.add_argument(
        "--includeSoftClipped",
        action="store_true",
        help=(
            "Include information from read bases that were marked as "
            "soft-clipped by the algorithm that made the BAM file. See "
            "https://samtools.github.io/hts-specs/SAMv1.pdf if this makes "
            "no sense to you."
        ),
    )

    parser.add_argument(
        "--compareWithPileupFile",
        metavar="FILE",
        help=(
            "Make a consensus using the pysam pileup function and compare "
            "it to the one made by usign pysam fetch. Write the summary to "
            'this file. If the file is "-", write to standard error.'
        ),
    )

    parser.add_argument(
        "--deletionSymbol",
        default="-",
        metavar="CHAR",
        help=(
            "When a deletion is detected, put this symbol into the "
            "consensus sequence. Use the empty string to have deletions "
            "omitted."
        ),
    )

    parser.add_argument(
        "--strategy",
        default="fetch",
        choices=("fetch",),
        help="The consensus-making strategy to use.",
    )

    parser.add_argument(
        "--deletionThreshold",
        default=0.5,
        type=float,
        metavar="N",
        help=(
            "If some reads have a deletion at a site and some do not, call "
            "the site as a deletion if the fraction of reads with the "
            "deletion is at least this value."
        ),
    )

    parser.add_argument(
        "--insertionCountThreshold",
        default=5,
        type=int,
        metavar="N",
        help=(
            "The number of reads that must have an insertion at an offset "
            "in order for the insertion to be called in the consensus."
        ),
    )

    parser.add_argument(
        "--progress",
        action="store_true",
        help=("Show a progress bar (unless standard error has been " "redirected)."),
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        help=(
            "Suppress diagnostic output. Note that this will silence "
            "warnings about differing reference names."
        ),
    )

    args = parser.parse_args()

    if args.version:
        print(__version__)
        sys.exit(0)

    try:
        consensus = consensusFromBAM(
            args.bamFilename,
            bamId=args.bamId,
            referenceFasta=args.referenceFasta,
            fastaId=args.fastaId,
            consensusId=args.consensusId,
            idLambda=args.idLambda,
            threshold=args.threshold,
            minCoverage=args.minCoverage,
            lowCoverage=args.lowCoverage,
            noCoverage=args.noCoverage,
            deletionSymbol=args.deletionSymbol,
            deletionThreshold=args.deletionThreshold,
            insertionCountThreshold=args.insertionCountThreshold,
            ignoreQuality=args.ignoreQuality,
            includeSoftClipped=args.includeSoftClipped,
            strategy=args.strategy,
            compareWithPileupFile=args.compareWithPileupFile,
            progress=args.progress,
            quiet=args.quiet,
        )
    except (ConsensusError, SamError) as e:
        print(f"{e} Exiting.", file=sys.stderr)
        sys.exit(1)
    else:
        print(consensus.toString("fasta"), end="")


if __name__ == "__main__":
    main()

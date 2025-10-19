#!/usr/bin/env python

import argparse
import os
import sys
from os.path import basename, join
from tempfile import mkdtemp

from dark.fasta import FastaReads
from dark.process import Executor

IVAR_FREQUENCY_THRESHOLD_DEFAULT = 0.6
IVAR_DOCS = "https://andersen-lab.github.io/ivar/html/manualpage.html#autotoc_md19"


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Make a consensus sequence.",
    )

    parser.add_argument("--reference", required=True, help="The reference FASTA file.")

    parser.add_argument(
        "--bam",
        help=(
            "The BAM file from which the consensus should be made. "
            "Required if --maskLowCoverage is used. If no BAM file is "
            "given, a VCF file must be provided. If both a BAM and a VCF "
            "file are given, the VCF file will take precedence."
        ),
    )

    parser.add_argument(
        "--vcfFile",
        help=(
            "The VCF file. If omitted, bcftools will be used to make a VCF "
            "file from the BAM file."
        ),
    )

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--id",
        help=(
            "The id to use in the consensus sequence in the output FASTA. "
            "If not given, the reference sequence id will be used."
        ),
    )

    group.add_argument(
        "--idLambda",
        metavar="LAMBDA-FUNCTION",
        help=(
            "A one-argument function taking and returning a read id. "
            "This can be used to set the id of the reference sequence based "
            "on the id of the reference sequence (the function will be "
            "called with the id of the reference sequence). E.g., "
            "--idLambda \"lambda id: id.split('_')[0]\" or "
            "--idLambda \"lambda id: id[:10] + '-consensus'\"."
        ),
    )

    parser.add_argument(
        "--sample",
        help=(
            "The name of the sample (from the @RG SM tag in the original "
            "alignment BAM file) for which a consensus should be made. "
            "If not given, the first sample name (from the #CHROM header) "
            "in the VCF file will be used."
        ),
    )

    parser.add_argument(
        "--dryRun",
        action="store_true",
        help="Do not run commands, just print what would be done.",
    )

    parser.add_argument(
        "--maskLowCoverage",
        default=0,
        type=int,
        help=(
            "Put an N into sites where the coverage is below the specified "
            "cutoff. If you specify a negative numer, masking will be "
            "turned off. Requires --bam."
        ),
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help=(
            "Show a log of commands that were (or would be, if --dryRun is "
            "used) executed."
        ),
    )

    parser.add_argument(
        "--noClean",
        action="store_false",
        dest="clean",
        help="Do not remove intermediate files or the temporary directory.",
    )

    parser.add_argument(
        "--callHaplotypesGATK",
        action="store_true",
        help=(
            "Use GATK to call haplotypes. See "
            "https://gatk.broadinstitute.org for details on GATK."
        ),
    )

    parser.add_argument(
        "--picardJar",
        help=(
            "The path to the Picard jar file. See "
            "https://github.com/broadinstitute/picard for details on "
            "Picard."
        ),
    )

    parser.add_argument(
        "--ivar",
        action="store_true",
        help="If given, ivar will be used to call the consensus.",
    )

    parser.add_argument(
        "--ivarFrequencyThreshold",
        type=float,
        help=(
            f"The frequency threshold used by ivar when calling the "
            f"consensus. If the frequency of the most-common nucleotide at "
            f"a site meets this threshold, the nucleotide will be called. "
            f"Otherwise, an ambiguous nucleotide code will be produced, "
            f"based on the smallest set of most-frequent nucleotides whose "
            f"summed frequencies meet the threshold. See {IVAR_DOCS} for "
            f"more information. If not given, "
            f"{IVAR_FREQUENCY_THRESHOLD_DEFAULT} is used. Can only be used "
            f"if --ivar is also specified."
        ),
    )

    parser.add_argument(
        "--ivarBedFile",
        help="If ivar should trim primers, a BED file of the primer positions.",
    )

    args = parser.parse_args()

    if not (args.bam or args.vcfFile):
        print("At least one of --bam or --vcfFile must be given.", file=sys.stderr)
        sys.exit(1)

    if args.maskLowCoverage and not args.bam:
        print("If --maskLowCoverage is used, --bam must be too.", file=sys.stderr)
        sys.exit(1)

    if args.ivar and not args.bam:
        print("If --ivar is used, --bam must be too.", file=sys.stderr)
        sys.exit(1)

    if args.ivarFrequencyThreshold is not None and not args.ivar:
        print(
            "If --ivarFrequencyThreshold is used, --ivar must be too.", file=sys.stderr
        )
        sys.exit(1)

    if args.ivar and args.ivarFrequencyThreshold is None:
        args.ivarFrequencyThreshold = IVAR_FREQUENCY_THRESHOLD_DEFAULT

    e = Executor(args.dryRun)

    tempdir = mkdtemp(prefix="consensus-")

    if args.vcfFile:
        vcfFile = args.vcfFile
    else:
        # No VCF file provided, so make one.
        vcfFile = join(tempdir, "vcf.gz")
        if args.callHaplotypesGATK:
            e.execute(f"samtools index {args.bam!r}")
            if args.picardJar:
                picardJar = args.picardJar
            else:
                try:
                    picardJar = os.environ["PICARD_JAR"]
                except KeyError:
                    print(
                        "If you use --callHaplotypesGATK, you must give a "
                        "Picard JAR file with --picardJar or else set "
                        "PICARD_JAR in your environment.",
                        file=sys.stderr,
                    )
                    sys.exit(1)

            indexFile = args.reference + ".fai"
            if os.path.exists(indexFile):
                removeIndex = False
            else:
                removeIndex = True
                e.execute(f"samtools faidx {args.reference!r}")

            if args.reference.lower().endswith(".fasta"):
                dictFile = args.reference[: -len(".fasta")] + ".dict"
            else:
                dictFile = args.reference + ".dict"

            if os.path.exists(dictFile):
                removeDict = False
            else:
                removeDict = True
                e.execute(
                    f"java -jar {picardJar!r} CreateSequenceDictionary "
                    f"R={args.reference!r} O={dictFile!r}"
                )

            e.execute(
                f"gatk --java-options -Xmx4g HaplotypeCaller "
                f"--reference {args.reference!r} "
                f"--input {args.bam!r} "
                f"--output {vcfFile!r} "
                f"--sample-ploidy 1 "
                f"-ERC GVCF"
            )

            if removeIndex:
                e.execute(f"rm {indexFile!r}")

            if removeDict:
                e.execute(f"rm {dictFile!r}")
        else:
            e.execute(
                f"bcftools mpileup --max-depth 5000 -Ou -f {args.reference!r} "
                f"{args.bam!r} | bcftools call --ploidy 1 -mv -Oz -o {vcfFile!r}"
            )

            e.execute(f"bcftools index {vcfFile!r}")

    if args.maskLowCoverage >= 0:
        # Make a BED file.
        bedFile = join(tempdir, "mask.bed")
        # The doubled-% below are so that Python doesn't try to fill in the
        # values and instead just generates a single % that awk sees.
        e.execute(
            f"samtools depth -a {args.bam!r} | "
            f"awk '$3 < {args.maskLowCoverage} "
            f'{{printf "%s\\t%d\\t%d\\n", $1, $2 - 1, $2}}\' > {bedFile!r}'
        )
        maskArg = "--mask " + bedFile
    else:
        maskArg = ""

    if args.sample:
        sample = args.sample
    else:
        result = e.execute(f"gunzip -c {vcfFile!r} | egrep -m 1 '^#CHROM' | cut -f10")
        sample = "SAMPLE-NAME" if args.dryRun else result.stdout.strip()

    consensusFile = join(tempdir, "consensus.fasta")

    if args.ivar:
        if args.ivarBedFile:
            tempBamFile = join(tempdir, basename(args.bam) + "-trimmed")
            result = e.execute(
                f"ivar trim -i {args.bam!r} -b {args.ivarBedFile!r} -p {tempBamFile!r} "
                "-e"
            )
            ivarTempBamFile = tempBamFile + ".bam"
            sortedIvarTempBamFile = tempBamFile + "-trimmed-sorted.bam"
            result = e.execute(
                f"samtools sort {ivarTempBamFile!r} -o {sortedIvarTempBamFile!r}"
            )
            bamFile = sortedIvarTempBamFile
        else:
            bamFile = args.bam

        ivarConsensusFile = join(tempdir, "temporary-consensus")
        result = e.execute(
            f"samtools mpileup -A -Q 0 {bamFile!r} | "
            f"ivar consensus -p {ivarConsensusFile!r} -q 20 "
            f"-t {args.ivarFrequencyThreshold!r} -m {args.maskLowCoverage!r}"
        )

        result = e.execute(f"mv {(ivarConsensusFile + '.fa')!r} {consensusFile!r}")

    else:
        result = e.execute(
            f"bcftools consensus --sample {sample!r} --iupac-codes {maskArg} "
            f"--fasta-ref {args.reference!r} {vcfFile!r} > {consensusFile!r}"
        )

        if not args.dryRun and result.stderr:
            print(result.stderr, end="", file=sys.stderr)

    if not args.dryRun:
        consensus = list(FastaReads(consensusFile))[0]
        if args.id is not None:
            consensus.id = args.id
        elif args.idLambda is not None:
            idLambda = eval(args.idLambda)
            consensus.id = idLambda(consensus.id)

        print(consensus.toString("fasta"), end="")

    if args.dryRun or args.log:
        print("\n".join(e.log), file=sys.stderr)

    if tempdir:
        if args.clean:
            e.execute(f"rm -r {tempdir!r}")
        else:
            print(f"Temporary directory {tempdir!r}.", file=sys.stderr)


if __name__ == "__main__":
    main()

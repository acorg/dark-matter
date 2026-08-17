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


def makeVcfWithGATK(args: argparse.Namespace, tempdir: str, e: Executor) -> str:
    vcfFile = join(tempdir, "vcf.gz")

    if args.picardJar:
        picardJar = args.picardJar
    else:
        try:
            picardJar = os.environ["PICARD_JAR"]
        except KeyError:
            sys.exit(
                "If you use --callHaplotypesGATK, you must give a Picard JAR file "
                "with --picardJar or else set PICARD_JAR in your environment."
            )

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

    return vcfFile


def makeVcfWithBcftools(args: argparse.Namespace, tempdir: str, e: Executor) -> str:
    vcfFile = join(tempdir, "vcf.gz")

    e.execute(
        f"bcftools mpileup --max-depth {args.maxDepth} -a AD,DP -q 20 -Q 20 -f "
        f"{args.reference!r} {args.bam!r} "
        f"| bcftools call --ploidy 1 -mv -Oz -o {vcfFile!r}"
    )

    e.execute(f"bcftools index {vcfFile!r}")

    filter_arg = "QUAL<20"

    if args.maskLowCoverage > 0:
        filter_arg += f" || FORMAT/DP<{args.maskLowCoverage}"

    filteredVcfFile = join(tempdir, "filtered-vcf.gz")
    e.execute(
        f"bcftools filter -e {filter_arg!r} {vcfFile!r} -Oz -o {filteredVcfFile!r}"
    )
    e.execute(f"mv {filteredVcfFile!r} {vcfFile!r}")
    e.execute(f"bcftools index -f {vcfFile!r}")

    return vcfFile


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
        "--maxDepth",
        default=5000,
        type=int,
        help=(
            "The maximum read depth to consider when calling the consensus via "
            "bcftools consensus."
        ),
    )

    parser.add_argument(
        "--maskLowCoverage",
        "--minimumCoverageDepth",
        "--minCoverageDepth",
        default=1,
        type=int,
        help=(
            "Put an N into sites where the coverage is below the specified "
            "cutoff. If you specify a negative numer, masking will be "
            "turned off. Note that specifying 0 is likely to cause problems, "
            "either with the bcftools or the samtools mpileup | ivar calling. "
            "If you're using bcftools and you set --maskLowCoverage to 0, "
            "bcftools will put the consensus base into all non-covered sites "
            "(because no low-coverage sites will be masked out by bcftools filter). "
            "If you're using ivar, the piped output from samtools mpileup -aa "
            "will include sites with zero coverage and the -m 0 given to ivar "
            "will result in ivar emitting nothing at all for those sites (so your "
            "consensus will end up being shorter than your reference, which may well "
            "not be what you want (i.e., might be totally wrong). Requires --bam."
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
        "--quiet",
        action="store_true",
        help="Do not emit a warning about --maskLowCoverage being set to zero.",
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

    if args.maskLowCoverage == 0 and not args.quiet:
        print(
            "WARNING: You have specified --maskLowCoverage 0. Do you for sure know "
            "what you're doing? If not, read the --help output. Use --quiet to "
            "suppress this warning.",
            file=sys.stderr,
        )

    if not (args.bam or args.vcfFile):
        sys.exit("At least one of --bam or --vcfFile must be given.")

    if args.maskLowCoverage and not args.bam:
        sys.exit("If --maskLowCoverage is used, --bam must be too.")

    if args.ivar and not args.bam:
        sys.exit("If --ivar is used, --bam must be too.")

    if args.ivarFrequencyThreshold is not None and not args.ivar:
        sys.exit("If --ivarFrequencyThreshold is used, --ivar must be too.")

    if args.ivar and args.ivarFrequencyThreshold is None:
        args.ivarFrequencyThreshold = IVAR_FREQUENCY_THRESHOLD_DEFAULT

    e = Executor(args.dryRun)

    tempdir = mkdtemp(prefix="consensus-")

    if args.bam:
        bamIndexFile = args.bam + ".bai"
        if os.path.exists(bamIndexFile):
            removeBamIndex = False
        else:
            removeBamIndex = True
            e.execute(f"samtools index {args.bam!r}")

    if args.vcfFile:
        vcfFile = args.vcfFile
    elif args.ivar:
        # We don't use a VCF file with iVar.
        vcfFile = None
    else:
        if args.callHaplotypesGATK:
            vcfFile = makeVcfWithGATK(args, tempdir, e)
        else:
            vcfFile = makeVcfWithBcftools(args, tempdir, e)

    consensusFile = join(tempdir, "consensus.fasta")

    if args.ivar:
        if args.ivarBedFile:
            tempBamFile = join(tempdir, basename(args.bam) + "-trimmed")
            e.execute(
                f"ivar trim -i {args.bam!r} -b {args.ivarBedFile!r} -p {tempBamFile!r} "
                "-e"
            )
            ivarTempBamFile = tempBamFile + ".bam"
            sortedIvarTempBamFile = tempBamFile + "-trimmed-sorted.bam"
            e.execute(f"samtools sort {ivarTempBamFile!r} -o {sortedIvarTempBamFile!r}")
            bamFile = sortedIvarTempBamFile
        else:
            bamFile = args.bam

        ivarConsensusPrefix = join(tempdir, "temporary-consensus")
        e.execute(
            # The samtools mpileup args were set on 2026-01-21 by Terry after looking at
            # the output of samtools mpileup (version 1.23).
            f"samtools mpileup -d 0 -aa -A -Q 0 {bamFile!r} | "
            f"ivar consensus -p {ivarConsensusPrefix!r} -q 20 "
            f"-t {args.ivarFrequencyThreshold!r} -m {args.maskLowCoverage!r}"
        )

        e.execute(f"mv {(ivarConsensusPrefix + '.fa')!r} {consensusFile!r}")

    else:
        sample = args.sample or (
            "SAMPLE-NAME"
            if args.dryRun
            else e.execute(f"bcftools query -l {vcfFile!r}").stdout.strip()
        )

        lowCoverage = join(tempdir, "low-coverage.bed")
        lowCoverageMerged = join(tempdir, "low-coverage-merged.bed")

        e.execute(
            f"samtools depth -a {args.bam!r} | "
            f"awk '$3<{args.maskLowCoverage} "
            f'{{printf "%s\\t%d\\t%d\\n", $1, $2 - 1, $2}}\' '
            f"| sort -k1,1 -k2,2n > {lowCoverage!r}"
        )
        e.execute(f"bedtools merge -i {lowCoverage!r} > {lowCoverageMerged!r}")

        e.execute(
            f"bcftools consensus --samples {sample!r} -m {lowCoverageMerged!r} "
            f"--fasta-ref {args.reference!r} {vcfFile!r} --output {consensusFile!r}"
        )

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

    if args.bam and removeBamIndex and args.clean:
        # We made the BAM index, so clean up.
        e.execute(f"rm {bamIndexFile!r}")


if __name__ == "__main__":
    main()

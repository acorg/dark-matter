#!/usr/bin/env python

import os
import sys
import argparse
import multiprocessing
from os.path import exists, join
from tempfile import mkdtemp
import pysam

from dark.process import Executor
from dark.bowtie2 import Bowtie2

DEFAULT_SAMTOOLS_VIEW_FLAGS = (
    pysam.FUNMAP | pysam.FSECONDARY | pysam.FDUP | pysam.FSUPPLEMENTARY
)


def saveStdin(args, e):
    if os.isatty(0):
        print("Reading sequences to match against from stdin.", file=sys.stderr)
    dirname = mkdtemp(prefix="run-bt2-stdin-", dir=args.tempdir)
    if args.tmpChmod:
        e.execute(f"chmod {args.tmpChmod} {dirname}")
    filename = join(dirname, "stdin.fastq")
    count = 0

    with open(filename, "w") as fp:
        for line in sys.stdin:
            fp.write(line)
            count += 1

    if count == 0:
        print("Standard input was empty! Exiting.", file=sys.stderr)
        sys.exit(1)

    return dirname, filename


def processMatch(args, e):
    """
    Run Bowtie2 to find matches.
    """
    if not args.index and not args.reference:
        print("One of --index or --reference must be given.", file=sys.stderr)
        sys.exit(0)

    if args.markDuplicatesPicard or args.callHaplotypesGATK:
        if args.picardJar:
            picardJar = args.picardJar
        else:
            try:
                picardJar = os.environ["PICARD_JAR"]
            except KeyError:
                print(
                    "If you use --markDuplicatesPicard or "
                    "--callHaplotypesGATK, you must give a Picard JAR file "
                    "with --picardJar or else set PICARD_JAR in your "
                    "environment.",
                    file=sys.stderr,
                )
                sys.exit(0)

    bt2 = Bowtie2(
        executor=e,
        threads=(multiprocessing.cpu_count() if args.threads is None else args.threads),
        verboseFp=(sys.stderr if args.verbose else None),
        dryRun=args.dryRun,
        tempdir=args.tempdir,
        tmpChmod=args.tmpChmod,
    )

    bt2.buildIndex(args.index or args.reference)

    if not args.align:
        if args.verbose:
            print("Bowtie2 alignment not done due to --noAlign.", file=sys.stderr)
        print("Bowtie2 temporary directory %r." % bt2.tempdir, file=sys.stderr)
        sys.exit(0)

    if args.fastq1:
        stdinDir, fastq1, fastq2 = None, args.fastq1, args.fastq2
    else:
        stdinDir, fastq1 = saveStdin(args, e)
        fastq2 = None

    if args.out:
        if exists(args.out) and not (args.force or args.dryRun):
            print(
                "Will not overwrite pre-existing output file %r. "
                "Use --force to make me." % args.out,
                file=sys.stderr,
            )
            sys.exit(1)

        if args.indexBAM:
            bai = args.out + ".bai"
            if exists(bai) and not (args.force or args.dryRun):
                print(
                    "Will not overwrite pre-existing output file %r. "
                    "Use --force to make me." % bai,
                    file=sys.stderr,
                )
                sys.exit(1)
            needBAI = True
        else:
            needBAI = False

    if args.callHaplotypesGATK:
        if args.vcfFile:
            for filename in (args.vcfFile, args.vcfFile + ".tbi"):
                if exists(filename) and not args.force:
                    print(
                        "Will not overwrite pre-existing VCF file %r. "
                        "Use --force to make me." % filename,
                        file=sys.stderr,
                    )
                    sys.exit(1)

    if args.callHaplotypesBcftools:
        if args.vcfFile:
            if exists(args.vcfFile) and not args.force:
                print(
                    "Will not overwrite pre-existing VCF file %r. "
                    "Use --force to make me." % args.vcfFile,
                    file=sys.stderr,
                )
                sys.exit(1)

    bt2.align(bowtie2Args=args.bowtie2Args, fastq1=fastq1, fastq2=fastq2)

    indexed = False

    # Sort, remove duplicates, etc. only if there is something in the BAM file.  Apart
    # from a slight efficiency gain, we check this because gatk dies with a typically
    # horrible Java stack dump if the are no reads in a BAM/SAM file.
    if not bt2.isEmpty():
        if args.sort and not (
            args.markDuplicatesPicard
            or args.markDuplicatesGATK
            or args.removePrimersFromBedFile
        ):
            bt2.sort()

        if args.removePrimersFromBedFile:
            bt2.sort()
            bt2.removePrimers(args.removePrimersFromBedFile)

        if args.markDuplicatesPicard:
            bt2.sort()
            bt2.picard(picardJar)

        if args.markDuplicatesGATK:
            bt2.sort(byName=True)
            bt2.markDuplicatesGATK()

        if args.removeDuplicates:
            bt2.removeDuplicates()

        if args.callHaplotypesGATK:
            bt2.makeBAM()
            bt2.indexBAM()
            indexed = True
            bt2.callHaplotypesGATK(
                picardJar=picardJar, vcfFile=args.vcfFile, referenceFasta=args.reference
            )

        if args.callHaplotypesBcftools:
            if not indexed:
                bt2.makeBAM()
                bt2.indexBAM()
                indexed = True
            bt2.callHaplotypesBcftools(vcfFile=args.vcfFile, referenceFasta=args.reference)

    if args.bam and args.indexBAM and not indexed:
        bt2.makeBAM()
        bt2.indexBAM()
        indexed = True

    if args.out:
        e.execute("mv '%s' '%s'" % (bt2.outputFile(), args.out))
        if needBAI:
            e.execute("mv '%s.bai' '%s.bai'" % (bt2.outputFile(), args.out))
    else:
        if not args.dryRun:
            with open(bt2.outputFile(), "rb") as fp:
                read, write = fp.read, sys.stdout.write
                while True:
                    data = read(2048)
                    if data:
                        write(data)
                    else:
                        break

    if stdinDir:
        # Remove the directory where we stashed standard input.
        e.execute("rm -r '%s'" % stdinDir)

    if args.clean:
        bt2.close()
    else:
        print("Bowtie2 temporary directory %r." % bt2.tempdir, file=sys.stderr)


def processOneIgnore(args, index, count, tempdir, e):
    """
    Process one ignored index.
    """
    if args.verbose:
        print("Preparing to ignore reads matching index %r." % index, file=sys.stderr)

    bt2 = Bowtie2(
        executor=e,
        threads=(multiprocessing.cpu_count() if args.threads is None else args.threads),
        verboseFp=(sys.stderr if args.verbose else None),
        dryRun=args.dryRun,
        tempdir=tempdir,
        tmpChmod=args.tmpChmod,
    )

    bt2.buildIndex(index)

    unAlignedFile = join(tempdir, "unaligned-%d" % count)

    if args.fastq1 and args.fastq2:
        bt2.align(
            bowtie2Args="%s --un-conc-gz %s" % (args.bowtie2Args, unAlignedFile),
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            discardSAM=True,
        )

        for i in "1", "2":
            src = unAlignedFile + "." + i
            dst = src + ".gz"
            e.execute("mv '%s' '%s'" % (src, dst))
            setattr(args, "fastq" + i, dst)

        if not e.dryRun:
            assert exists(args.fastq1)
            assert exists(args.fastq2)
    else:
        bt2.align(
            bowtie2Args="%s --un-gz %s" % (args.bowtie2Args, unAlignedFile),
            fastq1=args.fastq1,
            discardSAM=True,
        )

        e.execute("mv '%s' '%s.gz'" % (unAlignedFile, unAlignedFile))
        args.fastq1 = unAlignedFile + ".gz"
        if not e.dryRun:
            assert exists(args.fastq1)

    if args.clean:
        bt2.close()
    else:
        print(
            "Bowtie2 temporary ignore index directory %r." % bt2.tempdir,
            file=sys.stderr,
        )


def processIgnores(args, e):
    """
    Ignore the indices in args.ignoredIndices
    """
    if e.dryRun:
        tempdir = "/tmp/ignores"
    else:
        tempdir = mkdtemp(prefix="bt2-ignores-", dir=args.tempdir)
        if args.tmpChmod:
            e.execute(f"chmod {args.tmpChmod} {tempdir}")

    for count, ignoreIndex in enumerate(args.ignoredIndices):
        processOneIgnore(args, ignoreIndex, count, tempdir, e)
    return tempdir


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Run bowtie2 on a FASTA file. Optionally convert the "
            "result to BAM, sorting, and indexing."
        ),
    )

    parser.add_argument(
        "--index",
        help=(
            "Either: an accession number, a filename or the name of a "
            "pre-existing bowtie2 index (created with bowtie2-build). If "
            "not given and --reference is used, the reference will be "
            "used to build a bowtie2 index."
        ),
    )

    parser.add_argument(
        "--ignoreIndex",
        action="append",
        dest="ignoredIndices",
        help=(
            "Either: an accession number, a filename or the name of a "
            "pre-existing bowtie2 index (created with bowtie2-build). "
            "Reads matching this index will be ignored. May be repeated."
        ),
    )

    parser.add_argument(
        "--fastq1",
        "-1",
        help=(
            "The FASTQ reads to match against the bowtie2 index given by "
            "--index. Also use --fast2 if you have paired reads. "
            "If not given, single-end FASTQ reads will be read from "
            "standard input."
        ),
    )

    parser.add_argument(
        "--fastq2",
        "-2",
        help=(
            "The FASTQ reads to match against the bowtie2 index given by "
            "--index. Use this with --fastq1 to specify the mate "
            "file for paired-end reads."
        ),
    )

    parser.add_argument(
        "--bowtie2Args",
        default="--no-unal",
        help=(
            "Extra arguments to be passed to Bowtie2 (use --threads to "
            "specify a thread count)."
        ),
    )

    parser.add_argument(
        "--samtoolsViewArgs",
        default="-F %d -q 30" % DEFAULT_SAMTOOLS_VIEW_FLAGS,
        help="Arguments to be passed to samtools view to create the BAM file.",
    )

    parser.add_argument(
        "--tempdir",
        help=(
            "The temporary directory to use. If not specified, the value "
            "of the TMPDIR environment variable (if any) is used, or else "
            "/tmp."
        ),
    )

    parser.add_argument(
        "--out",
        "-o",
        help=(
            "The output file name. If not given, the resulting SAM or BAM "
            "will be written to standard output will be used."
        ),
    )

    parser.add_argument(
        "--reference",
        help=(
            "The reference FASTA file for use with --callHaplotypesGATK and "
            "--callHaplotypesBcftools. This will be used to build a Bowtie2 "
            "index if --index is not given."
        ),
    )

    parser.add_argument(
        "--vcfFile",
        help=(
            "The file to write VCF info to if --callHaplotypesGATK or "
            "--callHaplotypesBcftools are used."
        ),
    )

    parser.add_argument(
        "--markDuplicatesGATK",
        action="store_true",
        help=(
            "Use GATK to mark duplicates. See "
            "https://gatk.broadinstitute.org for details on GATK."
        ),
    )

    parser.add_argument(
        "--markDuplicatesPicard",
        action="store_true",
        help=(
            "Use Picard to mark duplicates. See "
            "https://github.com/broadinstitute/picard for details on "
            "Picard."
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
        "--tmpChmod",
        help=(
            "A chmod string for setting the permission on the temporary "
            "directory created. This will be passed to chmod(1), so you "
            'could specify "g+rwx", for example.'
        ),
    )

    parser.add_argument(
        "--removeDuplicates",
        action="store_true",
        help=(
            "Remove duplicates from the resulting SAM/BAM file. Best used "
            "in combination with an option that marks duplicates, such as "
            "--markDuplicatesGATK."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "Print a description of commands as they are (or would be, if "
            "--dryRun is used) executed."
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
        "--threads",
        type=int,
        help="The number of threads to use when running bowtie2 commands.",
    )

    parser.add_argument(
        "--noAlign",
        action="store_false",
        dest="align",
        help="Do not align with Bowtie2, just build an index.",
    )

    parser.add_argument(
        "--noBAM", action="store_false", dest="bam", help="Do not convert SAM to BAM."
    )

    parser.add_argument(
        "--noSort", action="store_false", dest="sort", help="Do not sort the BAM."
    )

    parser.add_argument(
        "--noIndexBAM",
        action="store_false",
        dest="indexBAM",
        help="Do not index the BAM file.",
    )

    parser.add_argument(
        "--noClean",
        action="store_false",
        dest="clean",
        help="Do not remove intermediate files or the temporary directory.",
    )

    parser.add_argument(
        "--force", action="store_true", help="Overwrite pre-existing output file."
    )

    parser.add_argument(
        "--dryRun",
        action="store_true",
        help="Do not run commands, just print what would be done.",
    )

    haplotypeCaller = parser.add_mutually_exclusive_group()

    haplotypeCaller.add_argument(
        "--callHaplotypesGATK",
        action="store_true",
        help=(
            "Use GATK to call haplotypes. See "
            "https://gatk.broadinstitute.org for details on GATK."
        ),
    )

    haplotypeCaller.add_argument(
        "--callHaplotypesBcftools",
        action="store_true",
        help="Use bcftools call to call haplotypes.",
    )

    parser.add_argument(
        "--removePrimersFromBedFile",
        help=(
            "If a bed file with Primers is specified, the Primers "
            "will be soft-clipped from the bam file using iVar"
        ),
    )

    args = parser.parse_args()

    if args.indexBAM and not args.bam:
        print(
            "The --indexBAM option only makes sense if you do not use " "--noBAM.",
            file=sys.stderr,
        )
        sys.exit(1)

    e = Executor(args.dryRun)

    if args.tempdir is None:
        args.tempdir = os.environ.get("TMPDIR", "/tmp")

    if args.ignoredIndices:
        ignoresDir = processIgnores(args, e)

    processMatch(args, e)

    if args.ignoredIndices:
        if args.clean:
            e.execute("rm -r '%s'" % ignoresDir)
        else:
            print(
                "Temporary directory with non-ignored inputs %r." % ignoresDir,
                file=sys.stderr,
            )

    if args.dryRun or args.log:
        print("\n".join(e.log), file=sys.stderr)


if __name__ == "__main__":
    main()

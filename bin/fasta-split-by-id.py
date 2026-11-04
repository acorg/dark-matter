#!/usr/bin/env python

import argparse
import sys
from math import log10
from os import mkdir, rename
from os.path import exists, join

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=(
        "Split sequences in a FASTA file into separate files, named "
        "either by their sequence id or numerically."
    ),
)

parser.add_argument("--outDir", default=".", help="The directory to make the files in.")

parser.add_argument(
    "--verbose",
    "-v",
    action="store_true",
    help="If given, print sequence ids as they are processed.",
)

parser.add_argument(
    "--numeric",
    action="store_true",
    help="If given, use numeric filenames, like 1.fasta, 2.fasta, etc.",
)

parser.add_argument(
    "--noLeadingZeroes",
    action="store_true",
    help="If given, numeric filenames will not have leading zeroes.",
)

parser.add_argument(
    "--filenamePrefix",
    default="",
    help="Use this as a prefix in output filenames.",
)

parser.add_argument(
    "--filenameSuffix",
    default="",
    help="Use this as a suffix (coming before '.fasta') in output filenames.",
)

parser.add_argument(
    "--idPrefix",
    default="",
    help="Use this as a prefix for sequence IDs.",
)

parser.add_argument(
    "--idSuffix",
    default="",
    help="Use this as a suffix for sequence IDs.",
)

parser.add_argument(
    "--force",
    "-f",
    action="store_true",
    help="If given, overwrite pre-existing files.",
)

parser.add_argument(
    "--dryRun",
    "-n",
    action="store_true",
    help="Don't create files, just print what would be done. Implies --verbose.",
)

parser.add_argument(
    "--saveAs",
    choices=("fasta", "fastq", "fasta-ss"),
    help=(
        "The output format. The default is to match the input format, "
        "so there is usually no need to specify this option. It can be "
        "used to force conversion from FASTQ to FASTA."
    ),
)

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)

if args.dryRun:
    args.verbose = True

if not exists(args.outDir):
    mkdir(args.outDir)

saveAs = (
    args.saveAs
    or (args.fasta and "fasta")
    or (args.fastq and "fastq")
    or (args.fasta_ss and "fasta-ss")
)

# Note: we may be reading the FASTA input from stdin, so we cannot read it
# more than once (and I don't want to store it all because it may be very
# large). That's why we do a second phase of processing to renumber the
# files we created if --numeric is used (and --noLeadingZeroes is not).

count = 0
for count, read in enumerate(parseFASTACommandLineOptions(args), start=1):
    idFields = read.id.split(maxsplit=1)
    id_ = idFields[0]
    base = (
        args.filenamePrefix
        + (str(count) if args.numeric else id_)
        + args.filenameSuffix
        + ".fasta"
    )
    filename = base if args.outDir == "." else join(args.outDir, base)
    read.id = (
        args.idPrefix
        + id_
        + args.idSuffix
        + (f" {idFields[1]}" if len(idFields) > 1 else "")
    )
    if args.dryRun:
        assert args.verbose
        if exists(filename) and not args.force:
            print(
                f"Would write {filename!r} ({len(read):,} nt) with ID "
                f"{read.id!r} but it already exists. Use --force to overwrite."
            )
        else:
            print(f"Would write {filename!r} ({len(read):,} nt) with ID {read.id!r}.")
    else:
        if exists(filename) and not args.force:
            print(
                f"Will not overwrite pre-existing file {filename!r}. "
                "Use --force to make me. Exiting.",
                file=sys.stderr,
            )
            continue
        if args.verbose:
            print(f"Writing {filename!r} ({len(read):,} nt) with ID {read.id!r}.")
        with open(filename, "w") as fp:
            print(read.toString(saveAs), end="", file=fp)


# Rename numeric filenames to have leading zeroes.
if count and args.numeric and not args.noLeadingZeroes:
    width = int(log10(count)) + 1
    if width > 1:
        if args.verbose:
            if args.dryRun:
                print("Renaming would be done to add leading zeroes")
            else:
                print("Renaming to add leading zeroes")
        for i in range(1, count + 1):
            old = str(i)
            new = "%0*d" % (width, i)
            if old == new:
                # We've reached the point where the new names have as many
                # digits as the old, so we can stop.
                break
            else:
                oldFilename = "%s.fasta" % join(args.outDir, old)
                newFilename = "%s.fasta" % join(args.outDir, new)
                if not args.force and exists(newFilename):
                    print(
                        f"Will not overwrite pre-existing file {newFilename!r}. "
                        "Use --force to make me. Exiting.",
                        file=sys.stderr,
                    )
                    continue

                if args.dryRun:
                    print(f"  {oldFilename!r} would be renamed {newFilename!r}")
                    continue

                if args.verbose:
                    print(f"  {oldFilename!r} -> {newFilename!r}")

                rename(oldFilename, newFilename)

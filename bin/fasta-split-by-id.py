#!/usr/bin/env python

import sys
import argparse
from os.path import join, exists
from os import mkdir, rename
from math import log10

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
    default=False,
    action="store_true",
    help="If given, print sequence ids as they are processed.",
)

parser.add_argument(
    "--numeric",
    default=False,
    action="store_true",
    help="If given, use numeric filenames, like 1.fasta, 2.fasta, etc.",
)

parser.add_argument(
    "--noLeadingZeroes",
    default=False,
    action="store_true",
    help="If given, numeric filenames will not have leading zeroes.",
)

parser.add_argument(
    "--force",
    default=False,
    action="store_true",
    help="If given, overwrite pre-existing files.",
)

parser.add_argument(
    "--saveAs",
    choices=("fasta", "fastq", "fasta-ss"),
    help=(
        "The output format. The default is to match the input format, "
        "so there is usually no need to specify this option. It can be "
        "used to force conversion from FASTQ to FASTA"
    ),
)

addFASTACommandLineOptions(parser)
args = parser.parse_args()
reads = parseFASTACommandLineOptions(args)

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
    id_ = read.id.split()[0]
    if args.numeric:
        base = "%d.fasta" % count
    else:
        id_ = read.id.split()[0]
        base = "%s.fasta" % id_
    if args.verbose:
        print("Writing", base)
    filename = join(args.outDir, base)
    if not args.force and exists(filename):
        print(
            'Will not overwrite pre-existing file "%s". Use --force to '
            "make me. Exiting." % filename,
            file=sys.stderr,
        )
        sys.exit(1)
    with open(filename, "w") as fp:
        print(read.toString(saveAs), end="", file=fp)

# Rename numeric filenames to have leading zeroes.
if count and args.numeric and not args.noLeadingZeroes:
    width = int(log10(count)) + 1
    if width > 1:
        if args.verbose:
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
                        'Will not overwrite pre-existing file "%s". Use '
                        "--force to make me. Exiting." % newFilename,
                        file=sys.stderr,
                    )
                    sys.exit(1)
                if args.verbose:
                    print("  ", oldFilename, "->", newFilename)
                rename(oldFilename, newFilename)

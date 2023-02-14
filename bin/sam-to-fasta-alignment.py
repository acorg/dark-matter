#!/usr/bin/env python

"""
Extract aligned (i.e., padded) queries in FASTA format from a SAM/BAM file.
"""

import sys
import argparse

from dark.filter import (
    addFASTAFilteringCommandLineOptions,
    parseFASTAFilteringCommandLineOptions,
)
from dark.reads import Reads
from dark.sam import SAMFilter, PaddedSAM
from dark.utils import nucleotidesToStr

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Produce aligned FASTA queries from a SAM/BAM file.",
)

parser.add_argument(
    "--rcSuffix",
    default="",
    help=(
        "A string to add to the end of query names that are reverse "
        "complemented. This is added before the /1, /2, etc., that are "
        "added for duplicated ids (if there are duplicates and "
        "--allowDuplicateIds is not used)"
    ),
)

parser.add_argument(
    "--rcNeeded",
    default=False,
    action="store_true",
    help=(
        "If given, queries that are flagged as matching when reverse "
        "complemented will be reverse complemented in the output. This "
        "must be used if the program that created the SAM/BAM input "
        "flags reversed matches but does not also store the reverse "
        "complemented query. The bwa program (mem and aln followed by "
        "samse) stores the queries reversed complemented if the match "
        "was, so this option is not needed for bwa. If in doubt, test the "
        "output of your matching program as this is very important!"
    ),
)

parser.add_argument(
    "--listReferenceInsertions",
    default=False,
    action="store_true",
    help=(
        "If given, information about reference sequence insertions will be "
        'printed to standard error. These correspond to "I" CIGAR '
        "operations that for the match would require inserting query bases "
        "into the reference. Because we cannot change the reference (in "
        "fact we typically do not have the reference in the SAM/BAM file), "
        "we cut the inserted bases out of the aligned query and save the "
        "information about what would have been inserted and where. That "
        "information is printed by this option. The output gives the "
        "0-based offset where the inserted base would be placed, followed "
        "by a list of the nucleotides that were suggested as being "
        "inserted and the number of times each nucleotide was suggested. "
        'So for example the output might contain "27: T:3, G:10" which '
        "indicates that 13 query (3 with T and 10 with G) matches would "
        "insert a nucleotide into the reference at offset 27."
    ),
)

SAMFilter.addFilteringOptions(parser)
addFASTAFilteringCommandLineOptions(parser)

args = parser.parse_args()
reads = parseFASTAFilteringCommandLineOptions(args, Reads())
samFilter = SAMFilter.parseFilteringOptions(args, filterRead=reads.filterRead)
paddedSAM = PaddedSAM(samFilter)

for read in paddedSAM.queries(rcSuffix=args.rcSuffix, rcNeeded=args.rcNeeded):
    print(read.toString("fasta"), end="")

if args.listReferenceInsertions:
    if paddedSAM.referenceInsertions:
        print(
            "(0-based) insertions into the reference:\n%s"
            % nucleotidesToStr(paddedSAM.referenceInsertions, "  "),
            file=sys.stderr,
        )
    else:
        print("No matches required an insertion into the reference.", file=sys.stderr)

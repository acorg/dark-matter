## 7.0.6 December 4, 2025

Added `--omitHeader`, `--ambiguous`, `--omitZeroes`, `--locations`, `--sort`,
`--python` to `bin/fasta-count-chars.py` to make it more flexible and
useful. Made it default to showing all DNA codes if `--chars` is not used (or
just the ambiguous codes if `--ambiguous` is given).

## 7.0.5 November 23, 2025

Change to make `newick-to-ascii.py` use the correct `ete4` method.

## 7.0.4 November 23, 2025

Upgrade `prseq` dependency to `0.0.33` to pick up its fix to fully read
multi-archive `gzip` and `bzip2` files.

## 7.0.3 November 13, 2025

Small fix to be compatible with `ete4`.

## 7.0.2 November 11, 2025

Un-skipped four tests due to fixing of https://github.com/bbuchfink/diamond/issues/905

## 7.0.1 November 11, 2025

Made run-bowtie2.py aware of 'large' bowtie index filenames (which end in
`.bt2l` not `.bt2`).

## 7.0.0 November 4, 2025

Upgrade from `ete3` dependency to `ete4`. The `curate-tree-ete3.py` script
has been renamed to `curate-tree-ete.py`. Given that the fix to upgrade is
so simple, I thought it better to break backward compatibility rather than
give the incorrect impression that `ete3` was still in use.

## 6.0.1 October 21, 2025

Make graphics reporting output go to stderr.

## 6.0.0 October 20, 2025

This release incorporates use of
[prseq](https://github.com/VirologyCharite/prseq/) for reading FASTA and
FASTQ files (see `FastaReads` in `src/dark/fasta.py` and `FastqReads` in
`src/dark/fastq.py`). The major version number has been incremented to 6
because this is a backwards-incompatible change, though the required code
adjustments will likely be very minor. The breaking change is that if you
were previously passing an (non-stdin) open file descriptor to `FastaReads`
or `FastqReads` it was from a file opened in text mode (from which Python
would read strings). In this release, if you want to pass an open file
descriptor to either of these classes, you must open your file in binary mode
(`mode="rb"`). That's because the open file descriptor will be passed to the
underlying Rust code (in `prseq`) and that needs to read bytes. If you are
not using `FastaReads` or `FastqReads` in this way, you are unlikely to
notice any change in this version. But... there are many `dark-matter`
scripts and functions that needed to be adjusted in minor ways in this
release, so it is possible that other tools that somehow depend on
`dark-matter` will also need to be modified in minor ways. Hopefully this
won't be too disruptive. I have absolutely no idea how many people actually
use this code. Please reach out on GitHub if you have questions/comments.

## 5.2.0 October 19, 2025

Make it possible to limit the number of reads shown in (CIV) blue plots. Save
FASTQ in gzipped format.

## 5.1.2 September 22, 2025

Added `--length` option to `bin/fasta-sort.py`.

## 5.1.1 September 19, 2025

Updated pyproject.toml so that src/dark subdirs are found.

## 5.1.0 September 6, 2025

Added `bin/msa-find-and-extract.py` to extract regions from multiple
sequence alignments (irrespective of gaps).  Removed `--includeShortGenomes`
and `--noShortGenomeWarning` from `bin/fasta-find.py`. Refactored
`bin/fasta-find.py` to be much cleaner. Added functions in `Read` and `Reads`
classes to find prefix/suffix delimited regions.

## 5.0.46 September 2, 2025

Added `--ignoreGaps`, `--gapCharacter`, and `--end` to `bin/fasta-find.py`.
The former makes it possible to search for sequences in gapped FASTA files
(e.g., those produced in a multiple-sequence alignment). The `--end` option
tells `bin/fasta-find.py` to report the offset of the end of the
match. Together these options make it possible to find the start and end of a
nucleotide sequence in gapped sequences.

## 5.0.45 August 26, 2025

Tiny bugfix to rotating a read for when the read sequence is empty.

## 5.0.44 August 18, 2025

Add gap canonicalization to `edlib` alignment output.

## 5.0.43 August 14, 2025

Added `mode` argument to `Reads.save` method to make it possible to append to
an existing file.

## 5.0.42 July 24, 2025

Slightly refactor `bin/compare-sequences.py`.

## 5.0.41 June 12, 2025

Third time lucky? Fix reporting of number of reads (as opposed to read ids)
in `bin/sam-coverage-depth.py`.

## 5.0.40 June 12, 2025

Fix recording of number of reads in `bin/sam-coverage-depth.py`.

## 5.0.39 June 10, 2025

Improve `bin/sam-coverage-depth.py` by a) being more flexible on finding the
reference sequence, b) ignoring (and warning about) non-ACGT nucleotides in
reads found by pysam in the BAM, and c) correcting the read counts.

## 5.0.38 June 4, 2025

Add printing of aa properties to `bin/aa-compare.py` (renamed from
`bin/codon-distance.py` which is now a symbolic link to `bin/aa-compare.py`).

## 5.0.37 May 10, 2025

Added `--includeDifferenceCounts` and `--includeDifferenceLocations` to
`bin/compare-sequences.py` and corresponding arguments to the underlying
function in `src/dark/dna.py` and tests in `test/test_dna.py`.

## 5.0.36 May 6, 2025

Move to using `pyproject.toml` and `uv`.

## 5.0.35 April 29, 2025

Added new `--minSiteMutationCount`, `--minBases`, and
`--maxIdenticalReadIdsPerSite` diversity filtering options to
`bin/sam-coverage-depth.py` and refactored the code therein to simplify
things. There should also be tests!

## 5.0.34 April 24, 2025

Generalized rerooting via `ete3` to allow multiple tip names to be given to
cause rerooting on the branch coming into the MRCA of those tips. Moved the
rerooting function into `dark/trees.py` and added tests.

## 5.0.33 April 21, 2025

Change name in `setup.py` to have an underscore following deprecation warning from PyPI.

## 5.0.32 April 21, 2025

Added `bin/remove-alrt-from-tree.py`.

## 5.0.31 March 10, 2025

Allow for zero reads matching a site when processing a pileup in `sam-coverage-depth.py`.

## 5.0.30 March 9, 2025

Add start/stop datetimes to summary output of `sam-coverage-depth.py`.

## 5.0.29 March 5, 2025

Added reference base and the 'from' base to the TSV output of `sam-coverage-depth.py`.

## 5.0.28 March 5, 2025

Added info to final output of `sam-coverage-depth.py`.

## 5.0.27 March 4, 2025

Many improvements to (and some simplifications of) `sam-coverage-depth.py`,
primarily a speed-up via parallelization.

## 5.0.26 Januery 27, 2025

Downgrade `dendropy` version requirement from `5.0.2` to `5.0.1`.

## 5.0.25 November 18, 2024

Adjust `bin/curate-tree-ete3.py` so it ensures the new root node is not formatted.

## 5.0.24 November 16, 2024

Made `bin/curate-tree-ete3.py` able to collapse clades due to low support or
alrt values, to remove support and/or alrt values, and to scale the length of
edges root leaving a new root.

## 5.0.23 November 15, 2024

Set the length on the new root node (when using `--detach` in calling
`bin/curate-tree-ete3.py`) to zero, seeing as that length is related to the
distance to the node that has been deleted (after rooting).

## 5.0.22 November 15, 2024

Added `--detach` option to `bin/curate-tree-ete3.py` to remove the outgroup
made by rooting on a node.

## 5.0.21 November 14, 2024

Added `bin/curate-tree-ete3.py`.

## 5.0.20 November 13, 2024

Small improvements to `bin/tree-info.py`.

## 5.0.19 November 12, 2024

Added simple `bin/tree-info.py` script to print information about tip names
and internal node labels and edge lengths in a phylogenetic tree.

## 5.0.18 November 9, 2024

Add explicit output format to samtools sort command in `dark/bowtie2.py`.

## 5.0.17 November 9, 2024

More messing with indexing BAM files.

## 5.0.16 November 9, 2024

Check if the SAM file in `dark/bowtie2.py` is empty by reading it and looking
for a non-header line. Undo the change of `5.0.15`. Added typing hints to
`dark/bowtie2.py`.

## 5.0.15 November 8, 2024

Always make BAM in `run-bowtie2.py`. Sigh.

## 5.0.14 November 8, 2024

Check whether a BAM/SAM file is empty in `run-bowtie2.py` in order to avoid
calling `gatk` to mark duplicates on a file with no mapped or unmapped reads,
since instead of just exiting gracefully, `gatk` crashes with the typical
Java runtime stack.

## 5.0.13 November 3, 2024

Added `bin/add-support-to-iqtree2-issue-343.py` script for adding support
labels to nodes in trees produced by `iqtree2` when (if not run with
`-keep-ident`) it adds nodes and tips for identical sequences. `iqtree2`
currently does not put a support label on (the edge leading to) the tip in
the original processing or onto the nodes introduced by adding tips for the
identical sequences. This is described in the `iqtree2` GitHub [issue
343](https://github.com/iqtree/iqtree2/issues/343).

## 5.0.12 November 1, 2024

This was merged late and became version `5.0.19`.

## 5.0.11 October 30, 2024

Added `--rotate` and `--maxWindows` options to `window-split-alignment.py`.

## 5.0.10 October 26, 2024

Added `intsToRanges` and `intsToStringRanges` to `dark/utils.py`.

## 5.0.9 October 24, 2024

Added `bin/curate-trees.py` script to collapse low-support branches in
phylogenetic trees (to make polytomies) and also to re-root and ladderize
them.

## 5.0.8 October 14, 2024

Added typing hints to `dark/process.py`.

## 5.0.7 October 5, 2024

Made `filter-fasta.py` print an error message when the `--checkResultCount`
check fails even if `--quiet` was used.

## 5.0.6 October 5, 2024

Add `rotate` method to `Read` class, plus tests.

## 5.0.5 October 5, 2024

Added `--rotate` option to `filter-fasta.py`.

## 5.0.4 October 5, 2024

Added `--upper`, `--lower`, `--upperId`, and `--lowerId` options to `filter-fasta.py`.

## 5.0.3 September 30, 2024

Removed `required` flag from `--out` option of `bin/plot-windowed-identity.py`.

## 5.0.2 September 28, 2024

Added `bin/plot-windowed-identity.py` along with `WindowedIdentity` class and tests.

## 5.0.1 September 22, 2024

Added lots of type hints. Completely removed six.

## 5.0.0 September 17, 2024

Removed `gb2seq` aligner options from `bin/sam-coverage-depth.py`. Bumped
major version number seeing as this will break things. I highly doubt anyone
is using the removed option. I (Terry) put it in for Christian Gabriel and
Annika Beyer during SARS-CoV-2 times.

## 4.0.89 August 7, 2024

Allow for `minWindow` to be `None` in `window-split-alignment.py`.

## 4.0.88 August 2, 2024

Added `bin/window-split-alignment.py` to split a FASTA alignment into windowed
sequences.

## 4.0.87 July 20, 2024

Added `stdout` and `stderr` options to `dark.process.Executor` to allow
for simple printing of executed commands and their outputs.

## 4.0.86 June 15, 2024

Added `bin/fasta-split-by-length.py`.

## 4.0.85 June 14, 2024

Added tiny `bin/reverse-complement.py` helper script.

## 4.0.84 June 1, 2024

Use a helper function (that puts the filename into failure exceptions) to
connect to sqlite3 databases.

## 4.0.83 May 23, 2024

Added `bin/ids.py` script to generate incrementing lists of ids.

## 4.0.82 May 23, 2024

Added `start` option to `dimensionalIterator`.

## 4.0.81 May 6, 2024

Added `asList` option to `parseRangeString` function in `utils.py`.

## 4.0.80 April 20, 2024

Added `--sequenceWhitelist`, `--sequenceBlacklist`,
`--sequenceWhitelistFile`, `--sequenceBlacklistFile`, `--sequenceRegex`,
and `--sequenceNegativeRegex` arguments to `filter-fasta.py` to match the
corresponding options for sequence titles (ids).

## 4.0.79 March 31, 2024

Added `--sortBy coverage` option to `bin/sam-reference-read-counts.py` and
moved much of its code into `dark/sam.py`.

## 4.0.78 March 15, 2024

Fix tiny bug in `bowtie2.py` so that BAM is not accidentally created when
removing duplicates.

## 4.0.77 March 15, 2024

Prevent an `IndexError` in `bin/sam-reference-read-counts.py` when there
are no matching reads.

## 4.0.76 March 14, 2024

Fixed (hopefully!) a bug in `bin/fasta-identity-table.py` that could cause
a `KeyError` when producing a non-square table,
[as described here](https://github.com/acorg/dark-matter/issues/780).

## 4.0.75 March 12, 2024

Moved `getNoCoverageCounts` from `bin/fasta-identity-table.py` into
`dark/reads.py`. Made it handle an empty set of no-coverage chars. Added
tests.

## 4.0.74 February 23, 2024

Added `allowedTaxonomicRanks` argument to the `SqliteIndexWriter` class for
building a protein database, and corresponding command line argument to
`bin/make-protein-database.py`.

## 4.0.73 January 20, 2024

Wrap the `savefig` call to make the pathogens panel (in
`dark/civ/proteins.py`) in a `try/except` to catch the `ValueError` that
results if the image has a dimension exceeding 2^16.

## 4.0.72 September 3, 2023

Made `aa-info.py` slightly more useful by identifying stop codons. Added
`--details` option to request printing of amino acid property numeric
details.

## 4.0.71 September 3, 2023

Corrected `aaVars` import in a few bins scripts.

## 4.0.70 August 30, 2023

Change to more recent `mysql-connector-python` in `setup.py`.

## 4.0.69 August 29, 2023

Added `bin/fastq-set-quality.py` script.

## 4.0.68 August 7, 2023

Fixed circular import.

## 4.0.67 August 7, 2023

Removed `bin/print-blast-xml.py`.

## 4.0.66 August 7, 2023

Removed `bin/print-blast-xml-for-derek.py`.

## 4.0.65 August 7, 2023

Nothing.

## 4.0.64 August 7, 2023

Removed `find-hits.py` from setup.py

## 4.0.63 August 7, 2023

Added type hints to scory.py and some to reads.py.

## 4.0.62 July 24, 2023

Removed `bin/find-hits.py`. Fixed matplotlib deprecation warnings. Ran black.

## 4.0.61 July 2, 2023

Backed out 4.0.59 and 4.0.60 changes.

## 4.0.60 July 2, 2023

More index checking in `bin/sam-coverage-depth.py`.

## 4.0.59 July 1, 2023

Subtract one from `stop` arg to pysam.pileup because, despite being
0-based, pysam apparently includes the final index.

## 4.0.58 June 12, 2023

Make `bin/download-refseq-viral-gbff.sh` exit non-zero if no viral genome
files are downloaded.

## 4.0.57 June 12, 2023

Use `[0-9]' instead of `\d` in `bin/download-refseq-viral-gbff.sh` so that
`egrep` works with the brew-installed GNU egrep on OS X.

## 4.0.56 June 8, 2023

Make `utils.py` `asHandle` work when passed a `pathlib.Path` argument.

## 4.0.55 May 6, 2023

Added `--includeAmbiguousMatches` and `--includeNonGapMismatches` options to
`bin/compare-sequences.py`.

## 4.0.54 May 6, 2023

Added `includeAmbiguousMatches` and `includeNonGapMismatches` options to
`matchToString` in `dark/dna.py`.

## 4.0.53 April 27, 2023

Replace use of `features` in `sam-coverage-depth.py` with the SAM file
reference length, to avoid trying to use an undefined `features` variable
when no reference is given.

## 4.0.52 April 24, 2023

Added `--force` option to `bcftools index` command in `dark/bowtie2.py`.

## 4.0.51 March 30, 2023

Added `--reverse` and `--complement` options to `bin/fasta-find.py` to tell
it to also look for simple inversions and complement sequences.

## 4.0.50 March 30, 2023

Added `bin/fasta-find.py` which, like `bin/fasta-match-offsets.py`, also
reports matching sequence offsets in FASTA files but can match numeric
regions and also look for reverse complemented matches.

## 4.0.49 February 19, 2023

Added `bin/fasta-match-offsets.py` to print offsets of sequence regular
expression matches.

## 4.0.48 February 14, 2023

More careful calling of `pysam.pileup` in `sam-coverage-depth.py` to avoid
an index error if a read mapping extends beyond the end of the reference
genome. `pysam` was returning an invalid `column.reference_pos` in that
case (invalid because the value is beyond the end of the reference, so it
can't be used as a reference offset).

## 4.0.47 February 10, 2023

Added printing of transversion and transition counts to
`bin/codon-distance.py`. Output is sorted first by distance, then by number
of transitions (highest to lowest) then by codon DNA. The idea being to
present the possible codon changes to get from one amino acid to another in
the order that requires the least change to the most.

The lists in the values in the dict returned by `codonInformation` now
contain 2-tuples instead of lists of length 2. If you pass
`countTransitions=True` to that function it will also put the count of the
number of transitions (as opposed to transversions) into the tuple. See the
tests in `test/test_codonDistance.py`.

## 4.0.46 January 16, 2023

Added `--reference` option to `sam-coverage-depth.py`.

## 4.0.45 November 27, 2022

Return an empty list of reads from mafft and needle when dryRun is True.

## 4.0.44 November 27, 2022

Added optional `executor` and `dryRun` arguments to MAFFT and needle
aligners to allow the caller to pass in a pre-existing process
executor. Added `--format` option to `newick-to-ascii.py` tree printing to
allow loading a wider range of Newick files.

## 4.0.43 Sept 25, 2022

Improve output of `compare-aa-sequences.py` to show the percentage of
matches in regions that do not involve a gap in either sequence.

## 4.0.42 Sept 21, 2022

Bump mysql connector version due to security issue with 8.0.13

## 4.0.41 Jul 9, 2022

Added `--sort` option to `bin/fasta-identity-table.py`.

## 4.0.40 Jul 5, 2022

Fixed small bug in `Reads.temporalBaseCounts`. Improved identity
calculation in `fasta-identity-table.py`.

## 4.0.39 Jul 2, 2022

Pass the `showNoCoverage` option value to the making of the HTML table in
`bin/fasta-identity-table.py`.

## 4.0.38 Jul 2, 2022

Added `--noNoCoverageLocations`, `--noCoverageChars`, and `--gapChars`
option to `bin/compare-sequences.py`. Fixed identity calculation bug in
`bin/fasta-identity-table.py` due to not including gaps resulting from the
pairwise alignment into the calculation.

## 4.0.37 Jul 2, 2022

Added `noCoverageChars` option to `compareDNAReads` and
`includeNoCoverageLocations` option to `matchToString`.

## 4.0.36 Jul 2, 2022

Make `compareDNAReads` more forgiving of unexpected characters in a DNA
sequence (specifically to deal with '?' that is used by Geneious to
indicate lack of coverage).

## 4.0.35 Jun 15, 2022

Added `--digits`, `--reverse`, `--sortBy`, `--header` and `--sortChars`
options to `bin/fasta-coverage.py`. Added `--regex` and `reverse` options
to `bin/fasta-sort.py`. Added `allowGaps` and `untranslatable` options to
`findORF` in `dark/reads.py`

## 4.0.34 Jun 15, 2022

Added simple `fasta-translate.py` script.

## 4.0.33 Jun 14, 2022

Small fix to text table output in `bin/fasta-identity-table.py`.

## 4.0.32 Jun 14, 2022

Added `--addZeroes` and `--highlightBest` options to
`bin/fasta-identity-table.py`.

## 4.0.31 Jun 14, 2022

Added `bin/fasta-coverage.py` command. Added `--upperOnly` option to
`bin/fasta-identity-table.py` (and fixed bug).

## 4.0.30 Jun 7, 2022

Added `temporalBaseCounts` method to `dark.Reads`.

## 4.0.29 June 3, 2022

Added `--align`, `--aligner`, and `--numberedColumns` options to
`bin/fasta-identity-table.py`. Added checking of pre-existing gap symbols
to the edlib aligner.

## 4.0.28 May 27, 2022

Added `-i` and `-v` options to `bin/genbank-grep.py`.

## 4.0.27 May 27, 2022

Renamed `bin/genbank-to-fasta.py` to `bin/genbank-grep.py` and made it able
to print in GenBank format as well as FASTA.

## 4.0.26 May 27, 2022

Added `bin/genbank-to-fasta.py` to extract FASTA from GenBank flat files.

## 4.0.25 May 27, 2022

Added `findORF` method to `DNARead` and `edlib` as an aligner option to
`compare-sequences.py`.

## 4.0.24 May 21, 2022

Added `--format` option to `ncbi-fetch-id.py` to allow fetching of GenBank
format flat files (use `--format gb`).

## 4.0.23 May 12, 2022

Removed a second unneeded viral check from `dark/taxonomy.py`.

## 4.0.22 May 12, 2022

Removed unneeded viral check from `dark/taxonomy.py`. Added 'generic' as a
pathogen type for protein reporting (in `dark/civ/proteins.py`) and that
can be given to `bin/proteins-to-pathogens-civ.py`.

## 4.0.21 April 13, 2022

Added `matchAmbiguous` option to `edlibAlign`.

## 4.0.20 April 4, 2022

Fixed `simple-consensus.py` argument error that shouldn't have been committed.

## 4.0.19 April 4, 2022

Added `edlib` alignment method.

## 4.0.18 February 9, 2022

Changed simple-consensus.py to use 'N' for no coverage and 'n' for low
coverage. Fixed issue with hard-clipped bases in CIGAR strings.

## 4.0.17 February 9, 2022

Made it so the genome/protein database builder can read compressed GenBank
or JSON files. Added simple scripts to download NCBI refseq viral FASTA or
GenBank flat file data.

## 4.0.16 February 5, 2022

Removed assertion of incorrect assumption that a read CIGAR string cannot
start with an insertion, based on a SAM file mapped using Geneious.

## 4.0.15 February 4, 2022

Added a simple consensus caller. Still needs to be stress tested in the
wild, and will need more tests and code refactoring to make it more
broadly useful.

## 4.0.14 January 20, 2022

Added internal optimization to make the SAM filtering fast when no
filtering is needed (this is the merge of a pull request from July 2021).

## 4.0.13 December 28, 2021

Added `bin/fasta-count-chars.py` script. Fixed argument bug in
`run-bowtie2.py` print statement.

## 4.0.12 November 29, 2021

Fixed another tiny bug when printing dry run output in
`bin/make-consensus.py`.

## 4.0.11 November 29, 2021

Fixed tiny bug when printing dry run output in `bin/make-consensus.py`.

## 4.0.10 November 2, 2021

Updated TravisCI config and README build status badge URL.

## 4.0.9 October 3, 2021

Moved some code out of dark/civ/proteins.py into dark/genbank.py to make it
more useful to others.

## 4.0.8 September 19, 2021

Added `DistanceMatrix` class to `dark/sam.py` for computing distances
between references according to which reads match both (and how well, if an
alignment score tag is given).

## 4.0.7 August 24, 2021

Color overall number of reads per pathogen in HTML output.

## 4.0.6 August 17, 2021

Fixed subtle bug introduced into `bin/filter-fasta.py` due to a code
reordering.

## 4.0.5 July 10, 2021

Added printing of reference lengths to `bin/sam-reference-read-counts.py`.
Changed the `--sort` option of that script to `--sortBy`. Added options to
stop printing references once they either have no reads mapped to them or
the number of new reads mapped to them ("new" as in not already mapping to
an earlier reference) falls to zero.

## 4.0.4 July 7, 2021

Added `fasta-split.py` command to split a FASTA/Q file into multiple files,
each containing a given number of sequences. Added `dark/utils.py` function
`take` to repeatedly yield lists of a given number of things from an
iterable.

## 4.0.3 June 15, 2021

Added `--topReferenceIdsFile` option to `sam-reference-read-counts.py` to
allow the ids of the best-matching reference to be saved. Probably this
should save the FASTQ.

## 4.0.2 June 15, 2021

Improved `sam-reference-read-counts.py` output to not double-count reads
that fall into multiple categories and also to report how many reads match
references that don't match any earlier reported reference (to give some
idea of how many reads uniquely match references, where 'unique' means
didn't match any other reported reference (with more reads when `--sort` is
also used).

## 4.0.1 June 14, 2021

Added `--sort` option to `sam-reference-read-counts.py` to sort output
(i.e., matched references) by highest number of mapped reads.

## 4.0.0 June 12, 2021

`btop2cigar` now returns a generator that gives the individual parts of a
CIGAR string. Use `''.join(btop2cigar(...))` to get the full CIGAR string.
This is not backwards compatible, hence the major version number change.

Pass the `tmpChmod` argument to a call to `Bowtie2` for which it was
missing.

## 3.2.7 April 3, 2021

Added `--maxMatchingReads` to title fitering.

## 3.2.6 April 3, 2021

Added BAM file sorting to `make-consensus.py` so ivar primer trimming actually
works.

## 3.2.5 March 4, 2021

Added ivar primer trimming (via the `--ivarBedFile` argument) to
`make-consensus.py`.

## 3.2.4 March 1, 2021

Changed how `ivarFrequencyThreshold` is set and checked in
`make-consensus.py`.

## 3.2.3 February 22, 2021

Added option to use `ivar` for consensus calling in `make-consensus.py`.

## 3.2.1 February 8, 2021

Added warning when `--removeDuplicatesUseMD5` is used without one of
`--removeDuplicates` or `--removeDuplicatesById`. Added `--tmpChmod` option
to `run-bowtie2.py`.

## 3.2.0 January 31, 2021

Dropped the `FastaFaiIndex` class. Use `SqliteIndex` instead. Fixed two
skipped Python 3 tests (for gzipped and bz2'd file reading).

## 3.1.88 February 02, 2021

Another attempt to fix the bug in bowtie2 that was introduced by adding the
`--removePrimersFromBedFile` option to `run-bowtie2.py`.

## 3.1.87 February 02, 2021

Fixed bug in bowtie2 that was introduced by adding the
`--removePrimersFromBedFile` option to `run-bowtie2.py`.

## 3.1.86 January 31, 2021

Fixed two failing tests due to changed DIAMOND (version 2.0.6) bitscore
calculation. If you have an older DIAMOND version you may need to update
it, if you plan to run the tests.

## 3.1.85 January 29, 2021

Added `--removePrimersFromBedFile` option to `run-bowtie2.py`. Fixed small
bug in `codon-distance.py`.

## 3.1.84 January 27, 2021

Added an option to trim primers using <b>ivar trim</b> to `run-bowtie2.py`.

## 3.1.83 August 26, 2020

`sam-coverage.py` now prints the min and max coverage per alignment too.

## 3.1.82 July 14, 2020

Added indexing to `callHaplotypesBcftools` in `bowtie2.py`.

## 3.1.81 July 13, 2020

Fixed errors from new version of flake8.

## 3.1.80 July 13, 2020

Added `--callHaplotypesBcftools` to `bowtie2.py`.

## 3.1.79 July 1, 2020

Added `--noFilter` option to `sam-coverage.py` and `sam-coverage-depth.py`
to allow them to run faster when no special filtering is needed.

## 3.1.78 May 25, 2020

Added `--callHaplotypesGATK` to `bin-make-consensys.py`.

## 3.1.77 May 22, 2020

Added code to combine multiple sequences (see `bin/combine-sequences.py`).

## 3.1.76 May 3, 2020

Changed `--maskNoCoverage` to `--maskLowCoverage` in `make-consensus.py` and
have it take an argument of the minimum coverage at which to call the
consensus.

## 3.1.75 April 27, 2020

Added `--sample-ploidy 1` to halpotype calling in `bowtie2.py`.

## 3.1.74 April 14, 2020

Improved `aa-info.py` to also match on partial full names. Improved
taxonomic detection of plant-only viruses.

## 3.1.73 March 27, 2020

Added `count` variable to `format-fasta.py` and a `--start` option to set
its starting value.

## 3.1.72 March 27, 2020

Added `format-fasta.py` script. Reverted unecessary hacks to
`fasta-sequences.py` to print MD5 sums, etc.

## 3.1.71 March 24, 2020

Added `--md5OneLine` argument to `fasta-sequences.py`. Prints TAB-separated
MD5 (sequence) sum, then the read id, then the sequence (and quality, if
any).

## 3.1.70 March 24, 2020

Added `--maxNFraction` argument to `filter-fasta.py`.

## 3.1.69 March 23, 2020

Added `--md5` arg to `fasta-sequences.py`.

## 3.1.68 March 16, 2020

Added expand all and collapse all buttons to HTML output generated by
`proteins-to-pathogens-civ.py`.

## 3.1.67 March 16, 2020

Tiny change to improve the HTML output generated by
`proteins-to-pathogens-civ.py`.

## 3.1.66 March 16, 2020

Many small changes to improve the HTML output generated by
`proteins-to-pathogens-civ.py`.

## 3.1.65 March 8, 2020

Fix idiotic logic error. Working too fast on complicated HTML-producing
code with no tests. But still idiotic.

## 3.1.64 March 8, 2020

Sort sample names in HTML output of `civ/proteins.py`. Added several names
to `PLANT_ONLY_VIRUS_REGEX` in `taxonomy.py`.

## 3.1.63 March 8, 2020

Added long comment and `ValueError` explanation regarding protein look-up
error likely due a protein database that's out-of-date with respect to
earlier result files.

## 3.1.62 March 8, 2020

Fix dumb error in accessing `proteinAccession` method.

## 3.1.61 March 8, 2020

Improve error messages to help debugging etc., when proteins cannot be
looked up and when exclusive host viruses are excluded in making a
protein/genome database.

## 3.1.60 March 8, 2020

Added `parseRangeExpression` to `utils.py`.

## 3.1.59 March 6, 2020

Fixed stupid argument mixup in calling `bcftools consensus` in
`make-consensus.py`.

## 3.1.58 March 6, 2020

Doubled some undoubled percentages in format string in `make-consensus.py`.

## 3.1.57 March 6, 2020

Added `--maskNoCoverage` option to `make-consensus.py`.

## 3.1.56 March 4, 2020

Added `sequenceToRegex` function to `dark/dna.py`.

## 3.1.55 March 3, 2020

Added `fasta-variable-sites.py` script.

## 3.1.54 Feb 21, 2020

Improve the printing of ambiguous sites to not show gaps in
`compare-sequences.py`.

## 3.1.53 Feb 21, 2020

Added `--showAmbiguous` to `compare-sequences.py`. Pass `--iupac-codes` to
`bcftools consensus` in `make-consensus.py`. Add `--max-depth 5000` to
`bcftools mpileup` call, also in `make-consensus.py`.

## 3.1.52 Feb 21, 2020

Make `sam-coverage-depth.py` not throw an error if there is no coverage at
all.

## 3.1.51 Feb 21, 2020

Pass the reference id to the `idLambda` function in `make-consensus.py`.

## 3.1.50 Feb 21, 2020

Moved `tempdir` assignment in `make-consensus.py` out one level so that it
happens including when we are given a VCF file.

## 3.1.49 Feb 21, 2020

Wrapped reads saving in a `Reads()` instance in `compare-sequences.py`.

## 3.1.48 Feb 21, 2020

Unreleased.

## 3.1.47 Feb 19, 2020

Added `--sampleName` and `--readGroup` options to `Bowtie2.align`. Added
`--id` and `--idLambda` options to `make-consensus.py` to make it possible
to set (or adjust) the name of the generated consensus.

## 3.1.46 Feb 19, 2020

Trivial change to setting of `tempdir` in `run-bowtie2.py`.

## 3.1.45 Feb 19, 2020

Change call to samtools.

## 3.1.44 Feb 19, 2020

Added `--tempdir` argument to `run-bowtie2.py`.

## 3.1.43 Feb 19, 2020

Added `make-consensus.py` script and tests (for a new `dark/bowtie2.py`
file). Upgraded `run-bowtie2.py`.

## 3.1.42 Feb 6, 2020

Added `--picard` to `run-bowtie2.py` command to allow removal of duplicates
found by Picard.

## 3.1.41 Feb 6, 2020

Added `run-bowtie2.py` command.

## 3.1.40 Jan 17, 2020

Make `compare-sequences.py` fall back to use `stretcher` if the call to
`needle` fails because the sequences are too long.

## 3.1.39 Jan 17, 2020

Fixed incorrect calculation of covered offset and total bases counts when
excluding reads based on minimum number of overlapping offsets in
`dark/genomes.py`.

## 3.1.38 Jan 13, 2020

Fix problem with bash set -u checking in `download-genbank.sh`.

## 3.1.37 Jan 13, 2020

Add final `genome-protein-summary.py` script.

## 3.1.35 Jan 13, 2020

Drop Python 3.8 from Travis and tox checking.

## 3.1.35 Jan 13, 2020

Drop Python2 from Travis and tox checking. Ugh.

## 3.1.34 Jan 13, 2020

Made sequence translation code work under Python 2 (again, even more
hopefully than the last time).

## 3.1.33 Jan 13, 2020

Made sequence translation code work under Python 2 (again, even more hopefully).

## 3.1.32 Jan 13, 2020

Made sequence translation code work under Python 2 (hopefully).

## 3.1.31 Jan 13, 2020

Added `sam-coverage-depth.py`.

## 3.1.30 Jan 13, 2020

Added `--minGenomeLength` to `make-protein-database.py`.

## 3.1.28 Dec 20, 2019

Removed the unused `taxonomy` (`VARCHAR`) column from the genomes table of
the protein database.

## 3.1.27 Dec 3, 2019

Fixed silly bug in alignment filtering code that was somehow not tested.

## 3.1.26 Dec 3, 2019

Added `--percentagePositiveCutoff` argument to `alignment-panel-civ.py`
and `noninteractive-alignment-panel.py` and all that implies, down to
reading DIAMOND output with the `ppos` value, storing it, restoring it,
and filtering on it in `dark.alignments.ReadsAlignmentsFilter`, plus tests.

## 3.1.25 Dec 3, 2019

Added `--minProteinCount` argument to `proteins-to-pathogens-civ.py` and
`proteins-to-pathogens.py`

## 3.1.24 Dec 3, 2019

Fixed a tiny Python <= 3.6 test output difference.

## 3.1.23 Dec 3, 2019

Fixed trivial Python2 incompatibility.

## 3.1.22 Dec 3, 2019

Added `--percentageIdenticalCutoff` argument to `alignment-panel-civ.py`
and `noninteractive-alignment-panel.py` and all that implies, down to
reading DIAMOND output with the `pident` value, storing it, restoring it,
and filtering on it in `dark.alignments.ReadsAlignmentsFilter`, plus tests.

## 3.1.21 Nov 30, 2019

Minor changes to HTML output.

## 3.1.20 Nov 30, 2019

Allow multiple `--preamble` args to `proteins-to-pathogens-civ.py`. Tiny
cosmetic adjustments to output HTML.

## 3.1.19 Nov 29, 2019

Added read count color levels indicator in HTML. Added
`--bootstrapTreeviewDir` to `proteins-to-pathogens-civ.py` and the `toHTML`
method of `dark.civ.proteins.ProteinGrouper`.

## 3.1.18 Nov 29, 2019

Fix (again) newline in HTML summary output.

## 3.1.17 Nov 29, 2019

Fix newline in HTML summary output.

## 3.1.16 Nov 29, 2019

Fix tiny bug in print arguments in `dark/civ/proteins.py`.

## 3.1.15 Nov 29, 2019

Added `--readCountColor` and `--defaultReadCountColor` to
`proteins-to-pathogens-civ.py` for differential coloring of read
counts. Added `citrus yellow` to plant-only virus name regex in
`dark/taxonomy.py`.

## 3.1.14 Nov 28, 2019

Added a const value for `--pathogenPanelFilename` in `proteins-to-pathogens-civ.py`.

## 3.1.13 Nov 26, 2019

Added `--dnaOnly` and `--maxGenomeLength` args to `make-protein-database.py`.
Improved `isRNAVirus` function so it returns `True` on retroviridae.
Improved the log output of the same script. Added a test (for HIV as an RNA
virus) and some small clean-ups.

## 3.1.12 Oct 21, 2019

Fixed incorrect `with` statement in `taxonomy.py`. Improved description of
fields in civ proteins HTML.

## 3.1.11 Sep 29, 2019

Standardized scripts that need a taxonomy database to use
`--taxonomyDatabase` on the command line and two utility functions in
`dark/taxonomy.py` to read them and also look in the
`DARK_MATTER_TAXONOMY_DATABASE` environment variable.

Added `dryRun`, `useStderr`, and handling of keyword arguments to the
`Executor.execute` method (in `dark/process.py`).

Add `LineageElement` to `taxonomy.py`. Get rid of `_preprocessLineage`
function and instead just have the `Taxonomy.lineageFromTaxid` method
adjust the 'no rank' ranks to be `-`. Added `skipFunc` and `stopFunc` to
lineage processing.

## 3.1.10 Sep 18, 2019

Make it so `get-taxonomy.py` and `get-hosts.py` can accept a name (e.g.,
Hymenoptera) as well as a taxonomy id or accession number.

## 3.1.9 Sep 12, 2019

Fixed `setup.py` error.

## 3.1.8 Sep 4, 2019

Added `get-hosts.py` and `get-taxonomy.py` scripts.

## 3.1.7 Sep 3, 2019

Added extremely basic `bin/describe-protein-database.py` command. To be
added to.

## 3.1.6 Sep 1, 2019

Added taxonomy info to HTML output.

## 3.1.5 Sep 1, 2019

Removed unguarded call to `self.pathogenPanel` that should have been
deleted on the last commit in `dark.civ.proteins.py` `toHTML` method.

## 3.1.4 Sep 1, 2019

Don't try to make a pathogen panel (in the `dark.civ.proteins.py` `toHTML`
method) if no pathogens were matched.

## 3.1.3 Sep 1, 2019

Fixed old code in `toStr`.

## 3.1.2 Sep 1, 2019

Fixed args in call to `toStr` in `proteins-to-pathogens-civ.py`.

## 3.1.1 Sep 1, 2019

Added `dark.civ` to packages in `setup.py`.

## 3.1.0 Aug 31, 2019

Added `make-protein-database.py`, `download-genbank.sh`, and
`parse-genbank-flat-file.py` scripts, as well as `doc/protein-database.md`
with some instructions on how to make a protein database. Added CIV
(Charite Institute of Virology) scripts `proteins-to-pathogens-civ.py` and
`alignment-panel-civ.py`.

## 3.0.80 Aug 19, 2019

Added `bin/create-newick-relabeling-output.py` and  `bin/relabel-newick-tree.py`.

## 3.0.79 Aug 8, 2019

Add `--omitVirusLinks` and `--omitSampleProteinCount` options to
`proteins-to-pathogens.py` to make HTML output less confusing when running
on RVDB or OKIAV databases.  Removed highlighting of pathogens with high
protein fraction since that was done in a non-useful way. Removed index
field from HTML output and removed HSP count unless it differs from the
read count.

## 3.0.78 Aug 2, 2019

Fixed silly import error.

## 3.0.77 Aug 2, 2019

Added link to per-pathogen reads in protein-to-pathogens HTML output for
Julia.

## 3.0.76 Jul 31, 2019

Slightly adjust appearance of HTML links for pathogens.

## 3.0.75 Jul 31, 2019

Added search link for ICTV for viruses in HTML output.

## 3.0.74 Jul 30, 2019

Added `--whitelistFile` and `--blacklistFile` options to
`noninteractive-alignment-panel.py`.

## 3.0.73 Jul 30, 2019

Adjust how protein and genome accession numbers are looked for in
`ProteinGrouper` depending on whether we guess they are NCBI or RVDB
sequence ids.

## 3.0.72 Jul 30, 2019

Make `NCBISequenceLinkURL` raise a more informative `IndexError` when it
cannot extract the wanted field.

## 3.0.71 Jul 30, 2019

Fixed stupid typo in `proteins-to-pathogens.py`.

## 3.0.70 Jul 30, 2019

Added `titleRegex` and `negativeTitleRegex` to `ProteinGrouper` class and
`--titleRegex` and `--negativeTitleRegex` arguments to
`proteins-to-pathogens.py`.

## 3.0.69 Jul 30, 2019

Added `--title` and `--preamble` args to output from
`proteins-to-pathogens.py`.  Fixed `ProteinGrouper` HTML NCBI protein link
and added genome link. Added positive and negative filtering by regex to
`TitlesAlignments` and tests. Improved NCBI link generation and tests.
>>>>>>> master

## 3.0.68 Jun 9, 2019

Refactored `SAMFilter` to allow filtering alignments in pileups. Added
`bin/sam-coverage.py`.

## 3.0.67 Apr 25, 2019

Use `dark.utils.StringIO` everywhere as it can be used as a context manager
in Python 2.

## 3.0.66 Apr 25, 2019

Added `pysam` to `setup.py` `install_requires` list. Removed `cffi`. Fixed
tests that were failing under Linux (apart from pyfaidx tests which are now
skipped on Linux). Removed mocking `File` class and replaced it with
`StringIO`.
>>>>>>> master

## 3.0.65 Jan 31, 2019

Fixed `AAread.ORFs` function in the `AARead` class and moved the
`--allowOpenORFs` (True/False) check to within the function. Added a
`DNAKozakRead` class. Changed `extract-ORFs.py` so that information
about Kozak consensus sequences can be returned.

## 3.0.64 Jan 28, 2019

Removed bone-headed use of full path to `fasta-join.sh` from
`bin/fasta-diff.sh`.

## 3.0.63 Jan 14, 2019

Added `compareAaReads` and `matchToString` to `aa.py`. Wrote tests in
`test_aa.py` for both. Moved `countPrint` to utils, used by `matchToString`
in `dna.py` and `aa.py`. Added `compare-aa-sequences.py` to the bin.

## 3.0.62 Dec 30, 2018

Added `matchToString` to `dna.py` to allow printing of a DNA match.

## 3.0.61 Dec 29, 2018

Added `--reverse` and `--reverseComplement` options to `filter-fasta.py`
and the underlying `ReadFilter` class.

## 3.0.60 Dec 13, 2018

In `reads.py`, changed the `_makeComplementTable` function so that
uppercase and lowercase bases are correctly reverse complemented into their
respective uppercase and lowercase complementary letters. Added a test to
`test/reads.py` to confirm that `reverseComplement` does this.

## 3.0.59 Dec 11, 2018

Added `--sampleName` option to `proteins-to-pathogens`.

## 3.0.58 Dec 6, 2018

Added `--maxORFLength` option to `extract-ORFs.py`. Fixed logic in
retrospect.

## 3.0.57 Dec 6, 2018

Added `--removeIds` option to `fasta-diff.sh` (and its helper script
`fasta-join.py`).

## 3.0.56 Dec 3, 2018

Make `convert-diamond-to-sam.py` print the correct (nucleotide) offset of
the start of the match, as though its subject sequence had been
nucleotides.

## 3.0.55 Dec 3, 2018

Make `convert-diamond-to-sam.py` print the list of required fields in case
an input line cannot be split into the expected number of fields.

## 3.0.54 Nov 27, 2018

`convert-diamond-to-sam.py` was putting the incorrect (AA, not nucleotide)
reference length into the SAM output. Introduced some spaces for easier
table layout into the HTML generated by `fasta-identity-table.py` when
called by `compare-consensuses.py`.

## 3.0.53 Nov 23, 2018

Fixed [#650](https://github.com/acorg/dark-matter/issues/650),
exception in `SAMFilter` when quality is `*` in a SAM file.

## 3.0.52 Nov 23, 2018

Added hard-clipping to CIGAR in SAM created by `convert-diamond-to-sam.py`.

## 3.0.51 Nov 23, 2018

Use `from six import StringIO` to avoid a PY2/3 incompatibility.

## 3.0.50 Nov 23, 2018

Added `bin/convert-diamond-to-sam.py` script to convert DIAMOND output
format 6 to SAM.

## 3.0.49 Nov 22, 2018

Added `btop2cigar` to `dark.btop` to convert BTOP strings to CIGAR strings.

## 3.0.48 Nov 12, 2018

Fixed `flake8` warnings about single backslashes in strings.

## 3.0.47 Nov 12, 2018

Make `fasta-diff.sh` use GNU parallel, if installed.

## 3.0.46 Nov 12, 2018

Make `fasta-diff.sh` handle compressed files, exit with diff's status,
and make it possible to pass command line options through to diff.

## 3.0.45 Nov 12, 2018

Added `bin/fasta-diff.sh` as a quick diff function that knows about
FASTA/FASTQ files.  Added `bin/fasta-join.py` as a helper function for
`bin/fasta-diff.sh`.

## 3.0.44 Nov 1, 2018

Fix [636](https://github.com/acorg/dark-matter/issues/636) in which SAM file
parsing threw an exception when an unmapped sequence with no CIGAR string
occurred in a SAM file (this can happen when running `bowtie2 --all`).

## 3.0.43 Nov 1, 2018

Fixed thinko in 3.0.42.

## 3.0.42 Nov 1, 2018

[Pysam issue 716](https://github.com/pysam-developers/pysam/issues/716)
wasn't solved the way we hoped it would be, so now the `filter-sam.py`
command always passes the template of the original SAM file to the
constructor for the new file. As a result, the new file will have `@SN`
entries for all the original sequences, even when `--referenceId` is passed
to `filter-sam.py` to restrict output to a specific set of reference ids.

## 3.0.41 Oct 31, 2018

Added `--showDiffs` option to `bin/compare-sequences.py`.

## 3.0.40 October 15, 2018

Force use of `mysql-connector-python` version `8.0.11` in
`requirements.txt` due to segmentation fault running tests using TravisCI.
Set version `0.5.0` of `pyfaidx` to avoid error in `pyfaidx/__init__.py`
line 711 (AttributeError: 'str' object has no attribute 'decode') in later
`pyfaidx` version. Removed some deprecation warnings when running tests.

## 3.0.39 October 10, 2018

The fix to solve [#630](https://github.com/acorg/dark-matter/issues/630)
was insufficient. That's fixed in this release, hopefully!

## 3.0.38 October 5, 2018

Fixed [#630](https://github.com/acorg/dark-matter/issues/630) to deal with
non-hard-clipped queries that have a CIGAR string that indicates they have
been clipped.

## 3.0.37 October 1, 2018

Add a `--titlesJSONFile` option to `noninteractive-alignment-panel.py`.

## 3.0.36 September 13, 2018

Added `whitelistFile` and `blacklistFile` to `ReadsAlignmentsFilter` class
in `dark/alignments.py`.

## 3.0.35 September 7, 2018

Fixed small bug in `filter-hits-to-fasta.py`.

Added flushing of intermediate output in `noninteractive-alignment-panel.py`.

## 3.0.34 July 18, 2018

Factored common SAM filtering code out into `dark.sam.SAMFilter`. Split
common FASTA command-line options into those for filtering (this is just
inclusion/exclusion) and those for editing (which change the FASTA records
in some way).

## 3.0.33 July 14, 2018

Added `compare-consensuses.py` script.

## 3.0.32 July 6, 2018

Added `bin/newick-to-ascii.py` script.

## 3.0.31 July 3, 2018

Added `storeQueryIds` option to `PaddedSAM.queries` method.

## 3.0.30 July 3, 2018

Added `alignmentCount` attribute to `PaddedSAM` class.

## 3.0.29 July 3, 2018

Renamed `alsoYieldAlignments` option of `PaddedSAM.queries` to `addAlignment`
and add the alignment to the `Read` instance instead of returning a tuple.

## 3.0.28 July 3, 2018

Added `alsoYieldAlignments` option to `PaddedSAM.queries` method to have
the returned generator also yield the `pysam.AlignedSegment` instance with
the gap-padded query sequence. This makes it possible to retrieve padded
queries from SAM/BAM and generate SAM/BAM (or FASTQ) of some subset of the
queries.

## 3.0.27 June 30, 2018

Added `bin/filter-sam.py` script.

## 3.0.26 June 26, 2018

Made a change in `dark/dna.py`, to make `identicalMatchCount` only count non-
ambiguous matches. Added testfunction `testMatchwithIdenticalAmbuguity`.

## 3.0.25 June 26, 2018

Added TravisCI Slack notifications.

## 3.0.24 June 25, 2018

Made `noninteractive-alignment-panel.py` option `--outputDir` to be required.
Added error message for this in `graphics.py`.

## 3.0.23 June 22, 2018

Changed the way reference sequence insertions are stored in a
`dark.sam.PaddedSAM` instance to make it possible to tell which query
sequences caused reference insertions.

## 3.0.22 June 21, 2018

Made `dark/sam.py` properly deal with secondary alignments that are missing
a SEQ.

## 3.0.21 June 21, 2018

Added `sam-reference-read-counts.py` script.

## 3.0.20 June 16, 2018

Updated ViralZone search URL in `dark/proteins.py`.

## 3.0.19 June 16, 2018

Added `--sites` argument to `compare-dna-sequences.py` and corresponding
`offsets` argument to the underlying function.

## 3.0.18 June 14, 2018

Fixed bug that got introduced when doing `3.0.17 June 14, 2018`.

## 3.0.17 June 14, 2018

Fixed bug that got introduced when doing `3.0.16 June 14, 2018`.

## 3.0.16 June 14, 2018

Made a change in `dark/proteins.py`, to make the `minProteinFraction` work
on a per sample basis, not per pathogen.

## 3.0.15 June 12, 2018

Fixed another bug (unreference variable) in `graphics.py` that crept in in
version `3.0.10 June 11, 2018`.

## 3.0.14 June 12, 2018

Fixed a bug in `diamond/alignments.py` that crept in in version
`3.0.10 June 11, 2018`.

## 3.0.13 June 11, 2018

Fixed a bug in `noninteractive-alignment-panel.py` that crept in after
version `3.0.10 June 11, 2018`.

## 3.0.12 June 11, 2018

* `pip install mysql-connector-python` now works, so added
`mysql-connector-python>=8.0.11` to `requirements.txt`, removed
`install-dependencies.sh`, and updated install line in `.travis.yml`.

## 3.0.11 June 11, 2018

* Added `bin/sam-to-fasta-alignment.py` script.

## 3.0.10 June 11, 2018

Dropped requirement that `noninteractive-alignment-panel.py` must be passed
information about the subject database. This is now only needed if
`--showOrfs` is given. The issue is that making the subject database can
take a long time and display of the subject ORFs is usuallly not needed.

## 3.0.9

Internal only.

## 3.0.8 June 1, 2018

* Added `--color` option to `fasta-identity-table.py`.

## 3.0.1 May 7, 2018

* Changed `Makefile` `upload` target rule.

## 3.0.0 May 5, 2018

* Moved all GOR4 amino acid structure prediction code into its own repo,
  at [https://github.com/acorg/gor4](https://github.com/acorg/gor4).
* As a result, the `gor4` method on the `dark.reads.AAread` class has been
  removed. This could be re-added by including `gor4` as a requirement but
  for now if you want that functionality you'll need to install `gor4`
  yourself and write a trivial function to call the `gor4` method on the
  read (or make a subclass of `AARead` that adds that method). I've done it
  this way because we have someone using the code who does not have a
  working C compiler and this was causing a problem building dark matter.
  Not a good reason, I know, but the GOR4 code makes for a good standalone
  code base in any case.

## 2.0.4 April 29, 2018

* Added `--sampleIndexFilename` and `--pathogenIndexFilename` to
  `proteins-to-pathogens.py`. These cause the writing of files containing
   lines with an integer index, a space, then a sample or pathogen name.
   These can be later used to identify the de-duplicated reads files for a
   given sample or pathogen name.

## 2.0.3 April 28, 2018

* Added number of identical and positive amino acid matches to BLAST and
  DIAMOND hsps.

## 2.0.2 April 23, 2018

* The protein grouper now de-duplicates read by id, not sequence.

## 2.0.1 April 23, 2018

* Fixed HTML tiny formatting error in `toHTML` method of `ProteinGrouper`
  in `dark/proteins.py`.

## 2.0.0 April 17, 2018

* The `--indices` option to `filter-fasta.py` was changed to accept a
  string range (like 10-20,25-30,51,60) instead of a list of single
  integers. It is renamed to `--keepSequences` and is also now 1-based not
  0-based, like its friends `--keepSites`.
* `--removeSequences` was added as an option to `filter-fasta.py`.
* The options `--keepIndices`, `--keepIndicesFile`, `--removeIndices`, and
  `removeIndicesFile` to `filter-fasta.py` are now named `--keepSites`,
  `--keepSitesFile`, `--removeSites`, and `removeSitesFile` though the old
  names are still supported for now.
* The `indicesMatching` methods of `Reads` is renamed to `sitesMatching`.
* `removeSequences` was added to read filtering in `ReadFilter` and as a
  `--removeSequences` option to `filter-fasta.py`.

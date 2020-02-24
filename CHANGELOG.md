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

Made sequence translation code work under Python 2 (again, even more hopefully than the last time).

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

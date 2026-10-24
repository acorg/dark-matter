# Categorizing the scripts in `bin`

Every script described below accepts a `--help` option, which prints a
full description of what it does along with its complete list of
command-line options. The one-line summaries below are necessarily
brief — if you're not sure exactly how a script behaves,
`script.py --help` is the place to check first.

## Running a script

**From a clone of this repository**, use [uv](https://docs.astral.sh/uv/)
from the repository's top-level directory:

```sh
uv run bin/fasta-lengths.py --help
```

This works the same way on Linux, macOS, and Windows.

**If you've installed dark-matter** (e.g., via `pip install dark-matter` or
`uv add dark-matter`, as described in the top-level `README.md`), the
scripts are installed into your Python environment's `bin` directory
alongside the library itself. If you're using uv, you don't need to find
that directory yourself — from within a project that depends on
dark-matter, `uv run` will find and run an installed script by name, on
any platform:

```sh
uv run fasta-lengths.py --help
```

Without uv, on Linux and macOS, that `bin` directory is normally on your
shell's `PATH` already (e.g., once you've activated a virtual environment
with `source .venv/bin/activate`), so you can just run a script by name:

```sh
fasta-lengths.py --help
```

On Windows, the installed scripts are plain `.py` files rather than `.exe`
launchers, so the reliable way to run one — regardless of your `PATH`
setup — is to invoke it explicitly with Python:

```bat
python path\to\your\venv\bin\fasta-lengths.py --help
```

(use `py` in place of `python` if that's what's set up on your system).

---

## 1. FASTA/FASTQ: general read inspection, filtering, and manipulation

* **`fasta-ids.py`** — Given FASTA on stdin, write the sequence ids to stdout.
* **`fasta-lengths.py`** — Given FASTA on stdin, write the sequence ids and lengths to stdout.
* **`fasta-count.py`** — Given FASTA on stdin, print the number of sequences to stdout.
* **`fasta-count-chars.py`** — Given FASTA on stdin, write information about the counts of characters in the sequences to stdout.
* **`fasta-join.py`** — Given FASTA on stdin, write the ids, sequences, and quality strings on a single TAB-separated line to stdout (used internally by `fasta-diff.sh`).
* **`fasta-sort.py`** — Write sorted FASTA/FASTQ to stdout. Sorting is by sequence id, then by sequence, then by quality (if FASTQ).
* **`fasta-split.py`** — Split sequences in a FASTA/FASTQ file into numerically named files, each containing a given number of sequences.
* **`fasta-split-by-id.py`** — Split sequences in a FASTA file into separate files, named either by their sequence id or numerically.
* **`fasta-split-by-length.py`** — Split FASTA/FASTQ sequences longer than a given length; print the beginning and end sequences of sufficiently long sequences separately.
* **`fasta-subset.py`** — Extract a subset of FASTA reads by id.
* **`fasta-subtraction.py`** — Given two or more FASTA file names as arguments, remove from the first file any read whose id also appears in a later file, and print what's left.
* **`filter-fasta.py`** — Given FASTA on stdin and a set of filtering criteria, write filtered FASTA to stdout. Has many options (see `--help`); the `bin/README.md` documents a couple of non-obvious uses (e.g., reformatting multi-line FASTA to one-line-per-sequence).
* **`filter-fasta-by-complexity.py`** — Filter FASTA for low-complexity sequences.
* **`filter-fasta-by-taxonomy.py`** — Filter FASTA based on taxonomy, matched against a taxonomic regular expression. **Note:** only works on FASTA with old-style NCBI `gi|...|` description lines, a numbering scheme NCBI retired in 2016 — likely a candidate for the "probably obsolete" pile (see bottom of this document).
* **`format-fasta.py`** — Given FASTA on stdin, write it to stdout according to a given output format.
* **`fasta-diff.py`** — Diff two FASTA sequences with visual highlighting (this is the Python-based, within-file diff; compare with the shell-based `fasta-diff.sh` below, which diffs whole FASTA/FASTQ *files*).
* **`fasta-diff.sh`** — A different, cruder FASTA/FASTQ diff: converts both files' records to single sorted TAB-separated lines and runs `diff`, then converts the result back into FASTA/FASTQ-like output. Good for a quick "are these basically the same" check, not for pinpointing exact differences.
* **`fasta-find.py`** — Given FASTA on stdin and a sequence (or genome region) to find, write information about (potentially overlapping) matches, optionally including reversed/complemented/reverse-complemented matches.
* **`fasta-match-offsets.py`** — Given FASTA on stdin and a sequence pattern, write information about which sequences the pattern matches, and where.
* **`fasta-base-indices.py`** — Given FASTA on stdin and a set of bases to look for, write the 1-based indices of where those bases occur in each sequence (or, with `--any`, where they occur in any sequence).
* **`fasta-identity-table.py`** — Print a FASTA sequence identity table (pairwise percent identity between all sequences).
* **`randomize-fasta.py`** — Read DNA FASTA from stdin and print FASTA to stdout, with the sequences randomized (ids and lengths preserved).
* **`fastq-set-quality.py`** — Given FASTA/FASTQ on stdin, write FASTQ to stdout with all quality scores set to a given value.
* **`combine-sequences.py`** — Combine multiple sequences into a single sequence.
* **`summarize-fasta-bases.py`** — Given FASTA on stdin, write a summary of sequence base categories to stdout (categories are currently hard-coded, not selectable on the command line).
* **`summarize-reads.py`** — Given a FASTA or FASTQ filename argument, print summary statistics: read count, total/average/median/min/max length, and per-base counts.
* **`fasta-coverage.py`** — Given a set of (presumably same-length or aligned) FASTA sequences, print per-character coverage/composition information (e.g., GC fraction) to stdout.
* **`trim-reads.py`** — Given a FASTA filename argument, detect a prefix/suffix common to all reads (e.g., leftover primer/adaptor) and print the reads with that trimmed off.
* **`trim-primers.py`** — Given FASTA on stdin, look for a primer sequence and write trimmed FASTA (the part after the primer) to stdout.
* **`split-fasta-by-adaptors.py`** — Given a set of known adaptors and FASTA on stdin, write separate FASTA files split by adaptor.
* **`adaptor-distances.py`** — Print a Levenshtein distance matrix for a set of known adaptors (helper for choosing/checking adaptors used by the two scripts above).
* **`reverse-complement.py`** — Print the reverse complement of one or more DNA sequences given directly as command-line arguments (not stdin).
* **`local-align.py`** — Given two FASTA files, perform a local Smith-Waterman alignment of all sequences in one against all those in the other.
* **`fasta-to-phylip.py`** — Given FASTA sequences of equal length on stdin, write Phylip format to stdout (for feeding into Phylip-family or other tree-building tools — could equally live in the "phylogenetic trees" category below).

## 2. DNA/amino-acid translation and codon tools

* **`dna-to-aa.py`** — Convert DNA to amino acids.
* **`aa-to-dna.py`** — Given amino acid FASTA on stdin, write FASTA with DNA sequences that would translate to the given amino acid sequences (codons chosen fixed or at random).
* **`fasta-translate.py`** — Translate FASTA and write the result to standard output.
* **`translate.py`** — Print the amino acid translation of one or more DNA sequences given as command-line arguments.
* **`extract-ORFs.py`** — Convert DNA or amino-acid FASTA to amino-acid, ORF-only FASTA.
* **`aa-info.py`** — Print detailed information (name, 1- and 3-letter abbreviations, codons, properties) for one or more amino acids or codons named on the command line, or for all amino acids if none are given.
* **`aa-compare.py`** — Compare two amino acids (or codons), given as command-line arguments: look each one up by name/abbreviation/codon, then print their shared and different properties and the possible codon-to-codon paths (with transition counts) between them, with colored terminal output.

## 3. Pairwise sequence comparison

* **`compare-sequences.py`** — Compare two DNA sequences: align them (MAFFT by default, or another aligner), then print a detailed match report (exact/ambiguous matches, mismatches, gaps, etc., via `matchToString`).
* **`compare-aa-sequences.py`** — The amino-acid equivalent of `compare-sequences.py`, aligning with `needle` (EMBOSS) and reporting the comparison via `dark.aa.matchToString`.

## 4. Multiple-sequence-alignment (MSA) and alignment-position analysis

* **`msa-find-and-extract.py`** — Find nucleotide sequences in an MSA and extract the matched region.
* **`msa-merge-gapped-region.py`** — Find nucleotide sequences bounding a region in an MSA and merge the gaps within that region into one gap, placed at a specified offset, moving the non-gap characters to one side (or around it).
* **`window-split-alignment.py`** — Split a FASTA alignment into equal-sized windows.
* **`plot-windowed-identity.py`** — Make a plot showing windowed nucleotide identity across an alignment.
* **`fasta-variable-sites.py`** — Given FASTA on stdin, write reads with only the variable sites to standard output (nothing is written if no sites are variable).
* **`position-summary.py`** — For a FASTA file of (presumably aligned) sequences, summarize what's happening at one specific (1-based) position.

## 5. Consensus-sequence building

* **`make-consensus.py`** — Make a consensus sequence from a FASTA reference and a BAM file, using external tools.

## 6. Read mapping and SAM/BAM analysis

* **`run-bwa.py`** — Run `bwa` on a FASTA file; optionally convert the result to BAM, sort, and index it.
* **`run-bowtie2.py`** — Run `bowtie2` on a FASTA file; optionally convert the result to BAM, sort, and index it.
* **`bwa-index-to-fasta.py`** — Extract FASTA from BWA index files (`.pac`/`.ann`/`.amb`).
* **`filter-sam.py`** — Given a SAM/BAM file and a set of filtering criteria, write filtered SAM/BAM to stdout.
* **`sam-merge.py`** — Concatenate, sort, and index SAM/BAM files with `samtools`.
* **`sam-references.py`** — Print SAM/BAM file reference names and lengths.
* **`sam-reference-read-counts.py`** — Print SAM/BAM file reference-sequence and unmapped-read counts.
* **`sam-coverage.py`** — Print SAM/BAM file coverage statistics.
* **`sam-coverage-depth.py`** — Print SAM/BAM per-site nucleotide composition and coverage statistics (sites with no coverage produce no output).
* **`plot-coverage-depth.py`** — Plot per-base coverage depth from a SAM/BAM file.
* **`sam-match-lengths.py`** — Print information about the lengths of query matches in a SAM/BAM file (query length, matched-region length excluding soft clips, soft-clip lengths on each end); can print per-query TSV and/or summary statistics.
* **`genome-protein-summary.py`** — Print SAM/BAM file protein-match statistics.
* **`extract-reference-from-fasta-and-bam.py`** — Extract a reference sequence from a FASTA file and its matches from a BAM file.
* **`reference-read-scores-to-JSON.py`** — Print a JSON object containing reference-to-read distances extracted from a SAM file.
* **`plot-references-by-inter-read-distance.py`** — Plot references according to their pairwise distances, computed from shared/non-shared read matches in a SAM file.

## 7. BLAST/DIAMOND hit processing and alignment panels

Scripts built around dark matter's BLAST/DIAMOND JSON record format (an
older, simplified JSON encoding of BLAST/DIAMOND hits used throughout this
part of the codebase).

* **`convert-blast-xml-to-json.py`** — Convert a BLAST XML file to dark matter's JSON format.
* **`convert-diamond-to-json.py`** — Convert DIAMOND tabular output to dark matter's JSON format.
* **`cat-json-blast-records.py`** — Read dark matter's JSON BLAST output (one JSON object per line) from stdin and pretty-print it to stdout.
* **`check-fasta-json-blast-consistency.py`** — Verify that a FASTA file and the JSON BLAST file produced from it are mutually consistent; exits non-zero with a stderr message if not.
* **`read-blast-json.py`** — Read dark matter's JSON BLAST records and report the elapsed time (a benchmarking/profiling script, not an end-user tool).
* **`read-blast-xml.py`** — Read dark matter's "simplified XML" BLAST records and report elapsed time (same purpose as above, for an XML-based variant).
* **`bit-score-to-e-value.py`** / **`e-value-to-bit-score.py`** — Convert between BLAST bit scores and e-values.
* **`filter-hits-to-fasta.py`** — Extract FASTA from matching BLAST hits.
* **`filter-reads-alignments.py`** — Given BLAST/DIAMOND JSON output and the FASTA/FASTQ files that were searched, print the *input reads* that pass a set of filtering criteria (the complement of `filter-hits-to-fasta.py`, which extracts matched *subjects* rather than queries).
* **`alignments-per-read.py`** — Read in BLAST records and print, for each read, a count of how many sequences it matches.
* **`get-reads.py`** — Print the reads that fall in a specified X-axis (read offset) and bit-score range of a BLAST alignment panel.
* **`alignment-panel-civ.py`** — Given BLAST/DIAMOND JSON output, the corresponding FASTA/FASTQ query files, and filtering criteria, produce a summary of matched titles and (optionally) an alignment panel image. The "civ" variants in this group are tailored to a specific subject-title format (`civ|GENBANK|<protein-id>|GENBANK|<pathogen-accession>|<product> [<pathogen name>]`) used by protein databases built with `make-protein-database.py`.
* **`proteins-to-pathogens-civ.py`** — Group protein match output (as produced by the alignment-panel scripts above) by the pathogen each protein comes from, using the pathogen accession number embedded in the "civ"-style subject title.

## 8. Protein/genome SQLite databases and GenBank data

* **`make-fasta-database.py`** — Create an sqlite3 database from FASTA sequences.
* **`make-protein-database.py`** — Create a genome/protein sqlite3 database from GenBank and JSON files; protein sequences for the database are printed to stdout for indexing by other tools. (This is the database consumed by the `-civ` scripts in category 7.)
* **`describe-protein-database.py`** — Print a description of a protein/genome database.
* **`subset-protein-database.py`** — Query an sqlite3 protein/genome database and write the matching sequences to stdout as FASTA.
* **`parse-genbank-flat-file.py`** — Parse GenBank flat files; exits non-zero if parsing fails (also used as a sanity check by `download-genbank.sh`).
* **`genbank-grep.py`** — Read a GenBank flat file and print GenBank or FASTA for matching records.
* **`download-genbank.sh`** — Download GenBank flat files for a list of accession numbers from NCBI, in parallel batches, with retries. One of a documented set of tools (see `doc/protein-database.md`) for building a protein database.
* **`download-refseq-viral-fasta.sh`** / **`download-refseq-viral-gbff.sh`** — Download the current RefSeq viral FASTA / GenBank flat-file release from NCBI's FTP site.
* **`ncbi-fetch-id.py`** — Fetch a sequence by id from NCBI and write it to stdout as FASTA.

## 9. Taxonomy

* **`get-taxonomy.py`** — Print taxonomy information (e.g., lineage) for given accession numbers or taxonomy ids, using a local taxonomy database.
* **`get-hosts.py`** — Print known hosts for given accession numbers or taxonomy ids, using a local taxonomy database.

## 10. Phylogenetic trees (Newick)

* **`curate-trees-dendropy.py`** — Curate phylogenetic trees (i.e., their specification/newick files) using `dendropy`.
* **`curate-tree-ete.py`** — The same kind of tree curation, but using the `ete4` library instead of `dendropy` (these look like two alternative implementations — worth deciding whether both are needed going forward, especially given the recent move to `ete4` noted in the CHANGELOG).
* **`newick-to-ascii.py`** — Print a Newick phylogenetic tree in ASCII.
* **`relabel-newick-tree.py`** — Relabel a Newick tree according to a file of from/to names.
* **`create-newick-relabeling-output.py`** — Find sequences with a specific property (via a user-supplied "relabel" function) and print a TAB-separated list of renamings, suitable as input to `relabel-newick-tree.py`.
* **`remove-alrt-from-tree.py`** — Read a Newick tree produced by `iqtree2` with both aLRT and bootstrap support values on its nodes, drop the aLRT values, optionally adjust the support values against given thresholds, and print the resulting tree.
* **`tree-info.py`** — Print information about phylogenetic tree nodes.

## 11. Small, general-purpose utilities

* **`ids.py`** — Print a generated sequence of ids (e.g., `id1`, `id2`, ... or a custom prefix and numeric range), with options for separators and zero-padding. General-purpose scripting/test-data helper; not specific to biological sequences.
* **`dark-matter-version.py`** — Print the installed `dark-matter` package version.

# Deprecated

* **`compare-consensuses.py`** — Compare different consensus-making methods (presumably running several of the above and diffing/scoring the results).
* **`filter-fasta-by-taxonomy.py`** - Depends on NCBI's GI numbering scheme, which NCBI retired in 2016, so these will silently do nothing useful on modern NCBI data unless you still have old-style FASTA/hit files lying around.
* **`filter-fasta-by-taxonomy.py`** — (also listed in category 1) filters FASTA by taxonomy, but only works with old GI-number description lines.
* **`noninteractive-alignment-panel.py`** — The same alignment-panel functionality as `alignment-panel-civ.py`, generated non-interactively (these two scripts look very close in structure/intent — worth checking whether one has simply superseded the other).
* **`proteins-to-pathogens.py`** — Group protein match output (as produced by the alignment-panel scripts above) by the pathogen each protein comes from, using the pathogen accession number embedded in the "civ"-style subject title.
* **`sam-to-fasta-alignment.py`** — Produce aligned FASTA queries from a SAM/BAM file.

# Incomplete

* **`simple-consensus.py`** — Make a consensus sequence directly from a BAM file using dark matter's own (pure-Python) consensus code, rather than external tools.
* **`convert-diamond-to-sam.py`** — Convert DIAMOND tabular output (produced with a specific required `--outfmt 6` field order, printable via `--printFields`) to SAM format.

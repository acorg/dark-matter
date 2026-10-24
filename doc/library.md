# The dark-matter library (`src/dark`)

This documents the functionality under `src/dark` for developers who
want to write their own code against the library, as opposed to
[`doc/scripts.md`](scripts.md), which documents the `bin` scripts for
people who just want to run existing tools. For each file: a short
overview, then a list of its public functions, classes (with their
important methods), and notable module-level data, plus a **Used by**
line naming which `bin` scripts *directly* import from that module
(not transitive usage through other `dark` modules).

---

## 1. Core reads & sequence I/O

The foundational `Read`/`Reads` class hierarchy and the
FASTA/FASTQ/PDB-FASTA readers built on it. See
[`doc/reads.md`](reads.md) for an explanation of the `Read` class
hierarchy (what each class is for and how they relate to one another)
— the summary below is deliberately brief by comparison.

### `src/dark/reads.py`

The core module of the library. Defines the `Read` class hierarchy (DNA, RNA,
amino-acid, secondary-structure, translated, and ORF reads), the `Reads`
collection class (and its filtering machinery), and the shared
`addFASTACommandLineOptions`/`parseFASTACommandLineOptions` helper pair used
by most `bin` scripts to read FASTA/FASTQ input.

- **`parseSiteSettings(settings, oneBased=False)`** — parses an iterable of `site:base[:quality]` strings (as used by the `--setSites` CLI option) into a list of `(offset, base, quality)` tuples, converting to 0-based offsets if `oneBased` is set. Raises `ValueError` on malformed input.
- **`Read`** — holds a single sequence: `id`, `sequence`, and an optional `quality` string (must match `sequence`'s length if given). Supports equality/ordering (`@total_ordering`), `len()`, hashing (based on an MD5 of id+sequence+quality), and slicing via `__getitem__` (returns a new instance of the same subclass).
  - `toString(format_="fasta")` / `toDict()` / `fromDict(d)` — convert to/from FASTA or FASTQ text, or a plain dict (for JSON-style serialization).
  - `setSites(settings, oneBased=False)` — mutates specific sequence (and quality) positions in place, given `site:base[:quality]` strings or pre-parsed tuples.
  - `reverse()` — reverses the sequence (not a complement).
  - `rotate(n, inPlace=True)` — rotates the sequence left/right by `n` characters.
  - `checkAlphabet(count=10)` — checks that the first `count` characters (or all, if `None`) are a subset of `self.ALPHABET`; raises `ValueError` if not.
  - `newFromSites(sites, exclude=False)` — returns a new read containing only (or, if `exclude`, excluding) the given 0-based site indices.
  - `find(pattern, start=0, end=False, caseSensitive=True, ignoreGaps=False, gapCharacter="-")` — finds a literal (non-regex) subsequence, optionally ignoring gap characters during the search while still returning an offset into the original (gapped) sequence.
  - `getPrefixAndSuffixOffsets(prefix, suffix, ...)` — finds the start offset of `prefix` and the end offset of `suffix` in the sequence (each -1 if not given/found); used by the prefix/suffix-trimming and region-extraction tooling.
  - `walkHSP(hsp, includeWhiskers=True)` — yields `(subjectOffset, residue, inMatch)` triples describing how this read aligns against a subject, as specified by an `HSP` (see `dark.hsp`); used for alignment-graph-style visualization.
  - `lowComplexityFraction()` — fraction of bases that are lowercase (the convention used elsewhere in the codebase for marking low-complexity regions).
- **`_NucleotideRead`** (internal base for `DNARead`/`RNARead`) — adds `reverseComplement()` and `translations()` (yields all 6 reading-frame `TranslatedRead`s).
- **`DNARead`** — `ALPHABET = {A,T,C,G}`; adds `findORF(offset, forward=True, requireStartCodon=True, allowGaps=True, untranslatable=None)`, which scans codon-by-codon from `offset` and returns a dict describing the ORF found (length, whether a start/stop codon was found, nucleotide sequence, and its translation).
- **`RNARead`** — `ALPHABET = {A,T,C,G,U}`; otherwise like `DNARead` but with an RNA complement table.
- **`DNAKozakRead`** — a `DNARead` subclass representing a detected Kozak consensus sequence; records the `originalRead`, `start`/`stop` offsets, and a `kozakQuality` percentage. Constructed by `findKozakConsensus` (in `dark.dna`).
- **`AARead`** — amino-acid reads. `properties()`/`propertyDetails()` map each residue to property flags/details (from `dark.aaVars`). `ORFs(openORFs=False)` scans the (already-translated) AA sequence for ORFs delimited by `M`...`*`, yielding `AAReadORF` instances; with `openORFs=True`, ORFs that run off either end of the sequence (no start or no stop codon found) are also yielded.
- **`AAReadWithX`** — `AARead` whose alphabet also allows `'X'`.
- **`AAReadORF`** — an `AARead` representing one ORF found within a larger `AARead`; tracks `start`/`stop` offsets and whether the ORF is `openLeft`/`openRight` (i.e., extends past the edge of the original read without a start/stop codon).
- **`SSAARead`** — an `AARead` paired with a `structure` string of the same length (PDB secondary-structure annotation); overrides `toString` to emit two FASTA records (sequence, then structure) in "fasta-ss" format.
- **`SSAAReadWithX`** — `SSAARead` combined with the `'X'`-extended alphabet.
- **`TranslatedRead`** — an `AARead` representing one DNA→AA translation of an `originalRead`, tagged with its `frame` (0/1/2) and whether it came from the `reverseComplemented` strand. `maximumORFLength(openORFs=True)` returns the longest ORF length found in it.
- **`ReadFilter`** — implements essentially all the filtering/editing logic behind `filter-fasta.py` and friends (length/complexity bounds, whitelist/blacklist by id or sequence (incl. regex and file-based lists), deduplication, random/fractional subsetting, site-based editing, id/sequence case changes, rotation, reverse/reverse-complement, quality-setting, arbitrary `idLambda`/`readLambda` functions, etc. — see its very long docstring for the full parameter list). Its `filter(read)` method applies all configured checks/edits in a fixed order and returns either the (possibly modified) read or `False` to drop it.
- **`readClassNameToClass`** — a `dict` mapping read class names (e.g. `"DNARead"`) to the classes themselves; used by `--readClass` CLI options and for deserialization. `getUnambiguousBases()` returns a similar mapping from class name to its set of unambiguous bases/residues.
- **`Reads`** — a lazily-iterable collection of `Read`s, combining: reads added via `add()`, an optional wrapped `initialReads` iterable, and reads yielded by the (overridable) `iter()` method (used by file-backed subclasses like `FastaReads`). Supports chained `.filter(**kwargs)` calls (each adds a `ReadFilter`), `save(filename, format_="fasta", mode="w")`, and numerous analysis methods:
  - `summarizePosition(index)` — per-site residue/base counts (`Counter`) across all reads.
  - `sitesMatching(targets, matchCase, any_)` — the set of sites where any (or all) reads have one of the `targets` bases.
  - `variableSites(confirm=False, homogeneityLevel=1.0, unknownAreAmbiguous=False)` — finds variable sites across equal-length reads, raising `ReadLengthsNotIdenticalError` otherwise.
  - `combineReads()` — collapses equal-length reads into one ambiguity-coded consensus string.
  - `temporalBaseCounts(firstPostId, ...)` — splits a time-ordered set of genomes into "pre"/"post" groups at a named id and reports newly-appearing bases and their frequencies in the "post" group.
  - `getPrefixAndSuffixOffsets(...)` / `getPrefixAndSuffixOffsetsForId(...)` / `extractRegion(...)` — locate a prefix/suffix across a set of (usually aligned) reads and extract the enclosed region (raising on conflicting offsets between reads).
  - `mergeGapsInRegion(...)` — within a prefix/suffix-delimited region, collapses all gap characters into one contiguous gap block (placed at `gapOffset`), used by `msa-merge-gapped-region.py`; can optionally write a per-sequence TSV summary and/or return a human-readable summary string.
  - `unfilteredLength()` — the true read count ignoring filters, but only available after the instance has been fully iterated once (raises `RuntimeError` otherwise).
- **`ReadsInRAM`** — a `Reads` subclass that eagerly consumes `initialReads` into a list, adding `len()`/indexing/assignment (`__len__`, `__getitem__`, `__setitem__`) on top of `Reads`' lazy iteration semantics.
- **`addFASTACommandLineOptions(parser)`** / **`parseFASTACommandLineOptions(args)`** — the shared argparse wiring used by ~30 scripts in `bin`. Adds a positional input-file argument (or `--fastaFile`, kept for backward compatibility; giving both is an error) plus `--readClass` and a mutually exclusive `--fasta`/`--fastq`/`--fasta-ss` group. If none of those three format flags is given, `_sniffFastaOrFastq` (private) peeks at the input's first non-whitespace byte (looking through gzip/bz2 compression if present) to auto-detect FASTA (`>`) vs FASTQ (`@`) without consuming the stream — raising `SystemExit` with a clear message if the format can't be determined this way (PDB FASTA, `--fasta-ss`, is never auto-detected since it's syntactically identical to plain FASTA). Returns the appropriate `FastaReads`/`FastqReads`/`SSFastaReads` instance.
- **`getNoCoverageCounts(reads, noCoverageChars)`** — counts, per read id, how many sequence characters are in the given no-coverage character set.
- **`simpleReadSplitter(length, leftPrefix="", leftSuffix="", rightPrefix="", rightSuffix="")`** — returns a function that splits reads longer than `length` into left/right fragments (with distinguishing id prefixes/suffixes), or passes shorter reads through unchanged; used by `fasta-split-by-length.py`.

**Used by:** `aa-to-dna.py`, `combine-sequences.py`, `compare-aa-sequences.py`, `compare-sequences.py`, `create-newick-relabeling-output.py`, `dna-to-aa.py`, `extract-ORFs.py`, `fasta-base-indices.py`, `fasta-count-chars.py`, `fasta-count.py`, `fasta-coverage.py`, `fasta-diff.py`, `fasta-find.py`, `fasta-identity-table.py`, `fasta-ids.py`, `fasta-join.py`, `fasta-lengths.py`, `fasta-match-offsets.py`, `fasta-sort.py`, `fasta-split-by-id.py`, `fasta-split-by-length.py`, `fasta-split.py`, `fasta-subset.py`, `fasta-to-phylip.py`, `fasta-translate.py`, `fasta-variable-sites.py`, `fastq-set-quality.py`, `filter-fasta.py`, `filter-hits-to-fasta.py`, `filter-sam.py`, `format-fasta.py`, `genome-protein-summary.py`, `msa-find-and-extract.py`, `msa-merge-gapped-region.py`, `plot-windowed-identity.py`, `randomize-fasta.py`, `reverse-complement.py`, `sam-coverage.py`, `sam-to-fasta-alignment.py`, `split-fasta-by-adaptors.py`, `summarize-fasta-bases.py`, `trim-primers.py`, `trim-reads.py`, `window-split-alignment.py`.

### `src/dark/fasta.py`

FASTA-specific reading, writing, de-duplication, and a disk-backed random-access index for FASTA sequences. Reading is done via the `prseq` package's `FastaReader`/`read_fasta`.

- **`fastaToList(fastaFilename)`** — Read an entire FASTA file (by filename) into a `list[Read]`, via `prseq.read_fasta`.
- **`dedupFasta(reads)`** — Given a `Reads` instance, yield only the first `Read` seen for each distinct sequence (duplicates are detected by an MD5 hash of the sequence string, not the id).
- **`dePrefixAndSuffixFasta(sequences)`** — Given an iterable of `Read`s, yield only those that are not a prefix or suffix (nor a duplicate) of some other, longer sequence in the input. Works by hashing every prefix and suffix of each sequence (longest sequences first) and skipping any sequence whose full hash has already been seen.
- **`fastaSubtract(fastaFiles)`** — Given a list of FASTA file names/handles, remove (by id) from the first file's reads any id that also appears in any of the subsequent files. Returns an iterator of the remaining `Read`s from the first file. Raises `IndexError` if given an empty list; asserts (don't just silently ignore) that a duplicate id has the same sequence in both files.
- **`FastaReads`** — A `Reads` subclass (see `dark/reads.py`) that reads one or more FASTA files/handles (or standard input, if `files=None`) and yields `Read` instances (of a configurable `readClass`, default `DNARead`) via its `iter()` method. `upperCase=True` uppercases each sequence as it's read.
- **`combineReads(filename, sequences, readClass, upperCase, idPrefix)`** — Build a single `Reads`/`FastaReads` by combining the contents of an (optional) FASTA file with a list of raw sequence strings given on the command line (each of which may embed its own id, by having the id as the text before the final space in the string; otherwise an id is auto-generated using `idPrefix` plus a counter). This is the standard way `bin` scripts that accept both a reference file and ad hoc sequences build their input.
- **`SqliteIndex`** — Build/use an sqlite3 database that maps FASTA sequence ids to `(file name, byte offset)` so that an arbitrarily large FASTA collection can be queried by id without loading it all into memory. (Unrelated to the differently-purposed `dark.civ.proteins.SqliteIndex`, which indexes genome/protein metadata rather than raw FASTA offsets.)
  - `addFile(filename)` — Scan (and index) all the sequences in a FASTA file (optionally BGZF-compressed, `.gz`/`.bgz`; plain `.bz2` is explicitly rejected since offset-seeking into a bzip2 stream isn't supported), recording each id's file number and byte offset. Raises `ValueError` on a duplicate file name or a duplicate sequence id (either within the file or across previously-added files).
  - `__getitem__(id_)` — Look up a sequence by id and read just that one sequence's bytes from disk by seeking to its recorded offset, returning it as a `Read` instance of the index's `readClass`. Raises `KeyError` for an unknown id.
  - `close()` — Close the underlying sqlite3 connection.

**Used by:** `alignments-per-read.py`, `alignment-panel-civ.py`, `check-fasta-json-blast-consistency.py`, `combine-sequences.py`, `compare-aa-sequences.py`, `fasta-subset.py`, `filter-fasta-by-taxonomy.py`, `fasta-diff.py`, `get-reads.py`, `local-align.py`, `make-consensus.py`, `make-fasta-database.py`, `filter-fasta-by-complexity.py`, `filter-reads-alignments.py`, `position-summary.py`, `filter-hits-to-fasta.py`, `randomize-fasta.py`, `noninteractive-alignment-panel.py`, `trim-primers.py`, `split-fasta-by-adaptors.py`, `sam-coverage-depth.py`.

### `src/dark/fastq.py`

A thin FASTQ-reading counterpart to `dark.fasta`.

- **`FastqReads`** — a `Reads` subclass. Constructor takes `files` (a single filename/handle, a list/tuple of them, or `None` for standard input) and a `readClass` (default `DNARead`) to instantiate for each record. `iter()` reads each file with `prseq.FastqReader` and yields `readClass(id, sequence, quality)` for every record.

**Used by:** `alignment-panel-civ.py`, `filter-reads-alignments.py`, `noninteractive-alignment-panel.py`.

### `src/dark/fasta_ss.py`

Reader for "PDB FASTA" files: pairs of FASTA records where the first holds an amino-acid sequence and the second (immediately following) holds its secondary-structure string of the same length.

- **`SSFastaReads`** — a `Reads` subclass. Constructor takes `_files` (filename/handle, list/tuple of them, or `None` for stdin — this `None`-handling was recently added so it matches `FastaReads`/`FastqReads`), a `readClass` (default `SSAARead`, must accept `(id, sequence, structure)`), and `upperCase` (uppercase both sequence and structure if set). `iter()` reads each file in sequence/structure record pairs via `prseq.FastaReader`, raising `ValueError` if a file has an odd number of records or if a sequence/structure pair's lengths differ, yielding `readClass(id, sequence, structure)` for each valid pair.

**Not used directly by any script in `bin`.**

### `src/dark/idutils.py`

A small generator-based utility for producing sequences of ids (e.g., `id1`, `id2`, ... or alphabetic/odometer-style ids), used to synthesize test or placeholder sequence ids.

- **`ids(start, end=None, prefix=None, sep=" ", zeroes=False, maxResults=None)`** — The main entry point: a generator that yields a run of ids from `start` up to (and including) `end`. `start`/`end` can each be a plain integer-like string (e.g., `"3"`) or an alphanumeric "odometer" string (e.g., `"aa"`, `"Az9"`), in which case each character position is incremented independently like an odometer wheel (digits cycle 0-9, lowercase a-z, uppercase A-Z; non-alphanumeric characters are held constant). If `end` is omitted, counting continues indefinitely from `start` unless `maxResults` caps it, or — if neither `end` nor `maxResults` is given — it is limited to incrementing only up to the next power of ten (to avoid an accidentally-infinite generator from plain `ids("3")`). If `prefix` is not given and both `start`/`end` are given, their longest common string prefix is factored out automatically and treated as a fixed prefix; `zeroes=True` zero-pads purely-numeric ids to the width of `end`.
- `_commonPrefix`, `_doInts`, `_doStrings` are private helpers implementing the above; not part of the public API.

**Used by:** `ids.py`.

### `src/dark/sequence.py`

Sequence-string utilities for locating a primer/adaptor within a sequence, forward and reverse-complemented.

- **`findPrimer(primer, seq)`** — Return a `list` of all (zero-based, non-overlapping) offsets in `seq` where `primer` occurs (case-insensitive).
- **`findPrimerBidi(primer, seq)`** — Return a `(forward, reverse)` tuple of two offset lists: where `primer` matches `seq` directly, and where the reverse complement of `primer` matches `seq` (i.e., where the original primer would match if `seq` were reverse complemented).
- **`findPrimerBidiLimits(primer, seq)`** — Given the forward/reverse offsets from `findPrimerBidi`, compute a `(start, end)` pair suitable for trimming all leading/trailing instances of the primer from a sequence, e.g. `seq[start:end]`: `start` is just past the last forward match (or `0` if none), and `end` is the first reverse match at or beyond `start` (or the first reverse match / full sequence length if there's no forward match).

**Used by:** `trim-primers.py`.

### `src/dark/summarize.py`

Two independent read-summarization helpers: overall per-file statistics, and a run-length-style categorization of a single sequence's characters.

- **`summarizeReads(file_handle, file_type)`** — Reads an entire FASTA (`file_type="fasta"`) or FASTQ (anything else) file via `FastaReads`/`FastqReads` and returns a `dict` with `read_number`, `total_length`, `average_length`, `max_length`, `min_length`, `median_length` (via `dark.utils.median`), and `base_counts` (a `defaultdict[str, int]` of base → total count across all reads). All numeric fields are `0` if there are no reads.
- **`sequenceCategoryLengths(read, categories, defaultCategory=None, suppressedCategory="...", minLength=1)`** — Walks a read's sequence, mapping each character through `categories` (falling back to `defaultCategory` for unmapped characters), and returns a `list` of `(category, runLength)` tuples for each maximal contiguous run of the same category. Runs shorter than `minLength` are not emitted individually; they're accumulated and emitted instead as a single `(suppressedCategory, totalSuppressedLength)` entry wherever they occur. Raises `ValueError` if `minLength < 1`.

**Used by:** `summarize-reads.py`, `summarize-fasta-bases.py`.

### `src/dark/simplify.py`

A single string-simplification helper for sequence titles/descriptions.

- **`simplifyTitle(title, target)`** — Given a title string and a `target` substring, return the prefix of `title` (as whole, space-separated words) up to and including the first word that contains `target` (case-insensitive); the matched word itself is truncated right after the end of the `target` match. If no word contains `target`, the whole title is returned unchanged. Intended for shortening verbose GenBank/NCBI-style sequence descriptions (e.g., trimming "Bovine polyomavirus DNA, complete genome" down to "Bovine polyomavirus" when `target="virus"`).

**Not used directly by any script in `bin`.**

---

## 2. DNA & amino-acid biology

Ambiguity-code handling, consensus-base calling, codon/amino-acid lookups and
properties, ORF-marker plotting, and SAM-CIGAR/BTOP alignment-encoding
conversions.

### `src/dark/dna.py`

Core DNA-sequence comparison, consensus, and ambiguity-code
logic. Defines the IUPAC ambiguous-nucleotide table and builds most of
the dark-matter DNA consensus/diff machinery on top of it.

- **`AMBIGUOUS`** — `dict[str, set[str]]` mapping each of the 15 IUPAC nucleotide codes (`A`, `C`, `G`, `T`, and the ambiguity codes `M`, `R`, `W`, `S`, `K`, `Y`, `V`, `H`, `D`, `B`, `N`) to the set of unambiguous bases it can represent. This is the central table nearly everything else in this module (and much of the rest of the DNA-handling code in the package) is built on.
- **`BASES_TO_AMBIGUOUS`** — The reverse mapping: a sorted, concatenated base-set string (e.g., `"ACGT"`) to its single ambiguity code (`"N"`).
- **`aaCompatible(s, aa)`** — Could the (possibly ambiguous) 3-character codon string `s` translate to amino acid `aa`? Checks both an exact `REVERSE_CODONS` hit and, failing that, whether `s`'s ambiguity is consistent with any of `aa`'s known codons.
- **`compatibleAAs(s)`** — Return the `set` of amino acids a (possibly ambiguous) codon `s` could translate to, by expanding all the concrete codons implied by `s`'s ambiguity codes and mapping each through `REVERSE_CODONS` (unknown/stop codons map to `"*"`).
- **`translatable(seq)`** — Is this nucleotide string cleanly translatable to protein (length a multiple of 3, no gaps/invalid bases, and no resulting `"X"` amino acid)? Returns `False` (never raises) for anything that fails.
- **`compareDNAReads(read1, read2, matchAmbiguous=True, gapChars="-", noCoverageChars=None, offsets=None)`** — The core pairwise DNA comparison function. Walks both sequences position by position (via `zip_longest`, so sequences of differing length are supported) and classifies each site as an identical match, an ambiguous match (if `matchAmbiguous`), a gap-involving mismatch, a "no coverage" site (characters in `noCoverageChars`, e.g. consensus sequences with uncovered regions), or a plain non-gap mismatch. `offsets`, if given, restricts comparison to just that set of positions. Returns a nested `dict` (see the function's docstring for the exact key layout) recording counts and offset-lists for all of the above, plus a `differenceOffsets` dict keyed by `"<char1><char2>"` pairs. This return value is the input to `matchToString`.
- **`matchToString(dnaMatch, read1, read2, matchAmbiguous=True, indent="", offsets=None, includeGapLocations=True, includeNoCoverageLocations=True, includeAmbiguousMatches=False, includeNonGapMismatches=False, includeDifferenceCounts=False, includeDifferenceLocations=False, includeNLocations=False)`** — Render the `dict` from `compareDNAReads` as a human-readable multi-line report string (exact/ambiguous match counts and percentages, mismatch breakdown by category, per-read gap/no-coverage/ambiguous/N counts and optional 1-based locations, and optional difference-count/location tables). The various `include*` flags control how much optional detail is appended. (This function's exact output format was recently revised — see `test/test_dna.py::TestMatchToString` for the current, authoritative expected-output examples.)
- **`findKozakConsensus(read)`** — Scan a `DNARead` for Kozak consensus sequences (the pattern `(gcc)gccRccATGG` around a start codon) and yield a `DNAKozakRead` (from `dark.reads`) for each `ATGG` start-codon context found, each annotated with a 0-100% "Kozak quality" score based on how many of the five variable positions match the preferred consensus base.
- **`FloatBaseCounts`** — Accumulate fractional base counts across a column of (possibly ambiguous) nucleotide codes, e.g. for building a consensus at one alignment position.
  - Constructor takes an iterable of codes; each ambiguous code's count is split fractionally across the bases it could represent (unknown codes count as fully-ambiguous `ACGT` if `unknownAreAmbiguous=True`, else as `"-"`).
  - `mostFrequent()` / `highestFrequency()` — the set of base(s) tied for the highest count, and that count's fraction of the total.
  - `homogeneous(level)` — is the most frequent base's fraction `>= level`?
  - `variable(confirm=True)` — is there more than one base present? With `confirm=True`, does a more careful check of whether the ambiguous codes could still all be consistent with a single real base (not just "multiple different symbols were seen").
  - `__len__`, `__str__` — number of distinct bases seen; a human-readable `"A:3 C:1 (0.750)"`-style summary.
- **`sequenceToRegex(sequence, wildcards="-?")`** — Convert a (possibly ambiguous) DNA sequence into a regular expression string, expanding each ambiguous code to a `[...]` character class; characters in `wildcards` that aren't otherwise ambiguity codes become `[ACGT]`. Raises `KeyError` on a genuinely unknown character.
- **`leastAmbiguous(bases)`** — Given an iterable of (unambiguous) base characters, return the single IUPAC code that exactly represents that set (via `BASES_TO_AMBIGUOUS`). Raises `KeyError` for an unknown base.
- **`leastAmbiguousFromCounts(bases, threshold)`** — Given a `dict`/`Counter` of base → frequency count and a homogeneity `threshold` (0-1), find the smallest set of (highest-count-first) bases whose cumulative frequency reaches `threshold`, including any bases tied with the last one included, then return the ambiguity code for that set via `leastAmbiguous`. This is the core of ambiguous-consensus-base calling. Raises `ValueError` for a negative count or threshold.
- **`Bases`** — A mutable accumulator of `(base, quality)` observations at one genome site, used when building a consensus from aligned reads.
  - `append(base, quality)` — record one observation (bases of `"N"` are ignored); returns `self` for chaining.
  - `consensus(threshold, minCoverage, lowCoverage, noCoverage)` — call the consensus base for this site: `noCoverage` if nothing was observed, `lowCoverage` if fewer than `minCoverage` observations were made, otherwise the ambiguous/unambiguous code from `leastAmbiguousFromCounts` at the given `threshold`.
  - Supports `+`/`+=` (merging two `Bases` instances), `__len__` (observation count), `__getitem__` (count for one base), and equality.

**Used by:** `compare-sequences.py`, `extract-ORFs.py`, `fasta-coverage.py`, `fasta-identity-table.py`, `fasta-count-chars.py`.

### `src/dark/aa.py`

The amino-acid analogue of `dna.py`: amino acid lookup/metadata and
pairwise amino-acid sequence comparison, built on the data tables in
`aaVars.py`.

- **`AminoAcid`** — A plain data holder for one amino acid's `name`, `abbrev3`, `abbrev1`, `codons` (tuple of its codons), `properties` (bitmask, see `aaVars.py`), `propertyDetails` (per-property scaled float values), and `propertyClusters` (per-property integer cluster numbers).
  - `propertySet()` — expand the `properties` bitmask into the `set` of individual property flag values it contains.
- **`find(s)`** — Look up amino acid(s) by a flexible specifier `s`: a full name (case-insensitively, with multi-word names handled specially), a 3- or 1-letter abbreviation, or a DNA/RNA codon (3 letters; `U`s are treated as `T`s). Falls back to substring-matching against full names if nothing else matches. Yields zero, one, or (for an ambiguous partial name match) more `AminoAcid` instances.
- **`propertiesForSequence(sequence, propertyNames, missingAAValue=-1.1)`** / **`clustersForSequence(sequence, propertyNames, missingAAValue=0)`** — For an `AARead`, extract, for each named property, a `list` of that property's (scaled float, or cluster int) value at every sequence position, substituting `missingAAValue` for unrecognized residues (e.g. `X`). Both delegate to the private `_propertiesOrClustersForSequence`, which validates that every requested property name is known (raising `ValueError` if not).
- **`matchToString(aaMatch, read1, read2, indent="", offsets=None)`** — The amino-acid counterpart of `dark.dna.matchToString`: render the `dict` returned by `compareAaReads` as a human-readable report (overall/non-gap match and mismatch counts and percentages, per-read gap counts and 1-based locations, extra-residue counts).
- **`compareAaReads(read1, read2, gapChars="-", offsets=None)`** — The amino-acid counterpart of `dark.dna.compareDNAReads`: a simpler site-by-site comparison (no ambiguity handling, since amino acid ambiguity codes aren't modeled here) producing match/mismatch/gap counts and per-read gap offsets, in the `dict` shape consumed by this module's `matchToString`.

**Used by:** `aa-compare.py`, `aa-info.py`, `compare-aa-sequences.py`.

### `src/dark/aaVars.py`

Pure data: the amino acid reference tables that `aa.py` (and, indirectly, much of the rest of the package) is built on. Essentially no logic beyond a few reverse-mapping dict/tuple comprehensions at import time.

- **`NAMES`** — 1-letter abbreviation → full amino acid name (e.g. `"A"` → `"Alanine"`). **`AA_LETTERS`** is the sorted list of its keys. **`NAMES_TO_ABBREV1`** is its reverse.
- **`ABBREV3`** — 1-letter → 3-letter abbreviation (e.g. `"A"` → `"Ala"`). **`ABBREV3_TO_ABBREV1`** is its reverse.
- **`HYDROPHOBIC`, `HYDROPHILIC`, `AROMATIC`, `SULPHUR`, `ALIPHATIC`, `HYDROXYLIC`, `TINY`, `SMALL`, `ACIDIC`, `BASIC_POSITIVE`, `NEGATIVE`, `POLAR`, `NONE`** — individual bit-flag constants (powers of two) for amino acid categorical properties. **`ALL_PROPERTIES`** is a tuple of all of them; **`PROPERTY_NAMES`** maps each flag to a human-readable name.
- **`PROPERTIES`** — 1-letter amino acid code → bitwise-OR of the categorical property flags above that apply to it (e.g. alanine is `HYDROPHOBIC | SMALL | TINY`).
- **`CODONS`** — 1-letter amino acid code → tuple of the (DNA) codons that encode it (per the standard genetic code). **`REVERSE_CODONS`** is the reverse mapping (codon → amino acid). **`START_CODON`** is `"ATG"`; **`STOP_CODONS`** is the 3-tuple of stop codons.
- **`PROPERTY_DETAILS`** — 1-letter code → `dict` of named physicochemical property → a *scaled* (-1 to 1) float value (aliphaticity, aromaticity, composition, hydrogenation, hydropathy, hydroxythiolation, iep, polar requirement, polarity, volume), drawn from several cited literature sources (see in-file comments for exact references). **`PROPERTY_DETAILS_RAW`** holds the same properties as their original, unscaled values.
- **`PROPERTY_CLUSTERS`** — 1-letter code → `dict` of the same named properties, but each mapped to a small integer cluster number (a coarser, discretized grouping of the raw/scaled values) rather than a continuous value.

**Used by:** `aa-compare.py`, `aa-to-dna.py`, `aa-info.py`, `summarize-fasta-bases.py`.

### `src/dark/codonDistance.py`

Small helper module for comparing codons: how many nucleotide
positions differ between two codons, and (for sets of codons encoding
two amino acids) the distribution of those distances, optionally
broken down into transitions vs. transversions.

- **`findDistance(co1, co2)`** — Returns the `int` number of positions (0–3) at which two length-3 codon strings differ (simple positional Hamming distance, no alignment).
- **`codonInformation(codons1, codons2, countTransitions=False)`** — Given two lists of codons (e.g., all codons for amino acid A and all codons for amino acid B), computes `findDistance` for every pair and groups the results. Returns a `dict` keyed by `int` distance (0–3), whose values are lists of `(codon1, codon2)` tuples, or `(codon1, codon2, transitionCount)` 3-tuples if `countTransitions=True` (a "transition" is a purine↔purine or pyrimidine↔pyrimidine substitution, i.e. one of A↔G or C↔T, per the module-level `TRANSITIONS` set).
- `TRANSITIONS` — module constant: a `set` of `(from, to)` nucleotide pairs considered transitions.

**Used by:** `aa-compare.py`.

### `src/dark/orfs.py`

Matplotlib plotting helpers that overlay start/stop codon markers onto
a sequence plot, split by reading frame. This is visualization support
code (not an ORF-finder in the sense of extracting ORF sequences — see
`dark.reads` for that).

- **`findCodons(seq, codons)`** — Generator over a `Bio.Seq.Seq` (or similar) scanning in consecutive, non-overlapping triplets from offset 0, yielding the `int` start offset of each triplet found in the given `codons` set. Note this only considers one fixed reading frame of `seq` itself; callers slice `seq[frame:]` first to scan the other two frames.
- **`addORFs(fig, seq, minX, maxX, offsetAdjuster)`** — Plots start (green) and stop (red) codon markers for all 3 forward reading frames of `seq` onto a matplotlib axes/figure `fig`, using `offsetAdjuster` (typically a `dark.intervals.OffsetAdjuster.adjustOffset`-like callable) to map sequence offsets to plot X coordinates. Sets axis limits/labels/title as a side effect.
- **`addReversedORFs(fig, seq, minX, maxX, offsetAdjuster)`** — Same as `addORFs` but for a reverse-complemented sequence, mirroring offsets via `maxX - offsetAdjuster(offset)`.
- `START_CODONS` / `STOP_CODONS` — module constants: `{"ATG"}` and `{"TAA", "TAG", "TGA"}` respectively.

**Used by:** Not used directly by any script in `bin` (it's `dark.graphics` that is the typical caller of these; no `bin` script imports `dark.orfs` directly).

### `src/dark/btop.py`

Parsing and conversion for BLAST/DIAMOND BTOP strings (a compact
alignment-difference encoding; format documented at
https://www.ncbi.nlm.nih.gov/books/NBK279682/). The main practical use
is turning a BTOP string into a SAM CIGAR string.

- **`parseBtop(btopString)`** — Generator that tokenizes a BTOP string into a stream of items: a bare `int` for a run of N identical/matching positions, or a `(queryLetter, subjectLetter)` 2-tuple for a single mismatch or gap (`"-"` denotes a gap in query or subject). Raises `ValueError` for malformed BTOP (e.g., a letter with no partner, two consecutive gaps, two identical consecutive letters, dangling trailing letter).
- **`countGaps(btopString)`** — Parses the string (via `parseBtop`) and returns a `(queryGaps, subjectGaps)` 2-tuple of `int` counts of `"-"` occurrences on each side.
- **`btop2cigar(btopString, concise=False, aa=False)`** — Generator yielding successive CIGAR-string pieces (e.g. `"12M"`) that together form a full CIGAR string (join with `"".join(...)`). Matches/mismatches become `=`/`X` (or both collapsed to `M` if `concise=True`); a query gap (`"-"` on the query side) becomes a reference deletion (`D`); a subject gap becomes an insertion (`I`). If `aa=True`, the BTOP numbers are interpreted as amino-acid counts and multiplied by 3 to produce nucleotide-length CIGAR runs (used for DIAMOND protein-vs-DNA alignments); in that mode you cannot distinguish an exact nucleotide match from a synonymous codon difference, so `concise` and `aa` cannot both be `True` (raises `ValueError` if so).

**Used by:** Not used directly by any script in `bin` (consumed via `dark.diamond.sam`, which calls `btop2cigar` internally).

### `src/dark/cigar.py`

Building and interrogating SAM CIGAR strings/operations, plus a couple
of helpers for working with pysam's `get_aligned_pairs`-style
(query-offset, reference-offset) pair lists. Imports
`dark.sam.CONSUMES_REFERENCE`, so this module sits "above" `dark.sam`
in the dependency order despite the name suggesting it's lower-level.

- **`dna2cigar(s1, s2, concise=False)`** — Builds a CIGAR `str` from two already-aligned, equal-length DNA strings (no gaps handled — this is purely match/mismatch, so it only ever emits `=`/`X` runs, or a single `M` run if `concise=True`). Raises `ValueError` if the sequences differ in length or are empty.
- **`makeCigar(reference, query, noEdgeInsertions=True)`** — Builds a CIGAR `str` from a pair of aligned strings that may contain `'-'` gap characters and leading/trailing space padding (space means "no base here for either side", used to represent one sequence extending beyond the other). Produces `M`/`I`/`D` runs; if `noEdgeInsertions=True` (the default), leading/trailing `I` runs are converted into soft-clips (`S`) instead, which is the convention SAM expects for a query that overhangs the reference. Raises `ValueError` if either input is empty.
- **`cigarTuplesToOperations(tuples, includeHardClip=True)`** — Generator that expands a pysam `cigartuples`-style list of `(operation, count)` pairs into a flat per-base sequence of operation codes (one yielded value per base), optionally eliding hard-clip (`CHARD_CLIP`) operations.
- **`softClippedOffset(offset, pairs, cigarOperations)`** — Given the index (`offset`) of a soft-clipped base within a pysam `get_aligned_pairs()` result (`pairs`) and the corresponding flat per-base `cigarOperations` list (e.g., from `cigarTuplesToOperations`), works out and returns the `int` reference coordinate that soft-clipped base would align to, by looking at the nearest non-hard-clipped neighbour (preferring looking backwards, then forwards). Uses `assert`s to validate its preconditions and raises `ValueError` if no usable neighbour is found.
- **`insertionOffset(offset, pairs, cigarOperations)`** — The insertion analogue of `softClippedOffset`: for a base that's part of an insertion (CIGAR `I`), finds the reference offset of the base immediately after the insertion. Returns a `(bool, int)` tuple — the `bool` says whether the answer was found by looking backward (and so should be read as "insertion occurs right after this reference position") vs. forward. Raises `ValueError` if no reference-consuming neighbour is found.
- Module constants `CINS_STR`, `CDEL_STR`, `CMATCH_STR`, `CEQUAL_STR`, `CDIFF_STR`, `CHARD_CLIP_STR` — the single-character SAM CIGAR operation codes (`I`, `D`, `M`, `=`, `X`, `H`) as strings, for use when building CIGAR text (as opposed to pysam's integer operation codes).

**Used by:** Not used directly by any script in `bin` (used internally by `dark.sam`, `dark.diamond.sam`, `dark.btop`, and `dark.local_align`).

### `src/dark/mutations.py`

Older, fairly ad hoc code for visualizing per-base mismatches of
BLAST-matched reads against a subject, and for studying APOBEC-style
mutation signatures (C→T / G→A hypermutation patterns with
trinucleotide context) in those mismatches. Several functions here
have hard-coded assumptions (e.g., about title string formats) that
make them specific to older BLAST-based workflows; this module
predates the DIAMOND/`dark.diamond` code path.

- **`basePlotter(blastHits, title)`** — For a given subject `title`, builds a textual "dot alignment" display: one line for the subject sequence, followed by one line per matching read where matching bases are replaced by `.` and mismatches show the actual (read) base. Requires `blastHits.plotParams`/`.titles[title]['plotInfo']` to already be populated (i.e. `computePlotInfo` must have been run on the `dark.blast.BlastHits`-like object first — asserts this). Returns a `list[str]`, each TAB-containing.
- **`getAPOBECFrequencies(dotAlignment, orig, new, pattern)`** — Scans a `basePlotter`-style dot alignment for positions where the subject had base `orig` and a read had `new` there, and tallies those occurrences by the surrounding trinucleotide context, restricted to one of two fixed context sets (`pattern` must be `'cPattern'` or `'tPattern'`, i.e. only contexts with `C` or `T` in the middle, matching APOBEC3 substrate motifs). Returns a `dict` (really a `Counter`-like `defaultdict(int)`) mapping trinucleotide context → count.
- **`getCompleteFreqs(blastHits)`** — Drives `basePlotter` + `getAPOBECFrequencies` across every title in `blastHits` and all 6 possible single-base substitution types (`C>A`, `C>G`, `C>T`, `T>A`, `T>C`, `T>G`), returning a nested `dict`: `{title: {substitutionType: {context: count}, "numberOfReads": int, "bitScoreMax": float}}`.
- **`makeFrequencyGraph(allFreqs, title, substitution, pattern, color='blue', createFigure=True, showFigure=True, readsAx=False)`** — Draws (via matplotlib) a bar chart of one substitution type's context frequencies (normalized by read count) for one title, either as a standalone figure or into a supplied subplot `readsAx`. Returns the `float` max bar height used, for axis-scaling by callers.
- **`makeFrequencyPanel(allFreqs, patientName)`** — Builds a grid of `makeFrequencyGraph` subplots (one row per title, one column per substitution type) into a single matplotlib figure, normalizing all Y axes to a common scale, and returns `allFreqs` unchanged (the useful side effect is the displayed figure).
- **`mutateString(original, n, replacements='acgt')`** — Returns a new `str` with exactly `n` randomly-chosen positions of `original` replaced by a randomly-chosen different character from `replacements`. Raises `ValueError` if `original` is empty, `n` exceeds `len(original)`, `replacements` has duplicate characters, or `replacements` has length 1 and already appears in `original` (making replacement impossible at some position).

**Used by:** Not used directly by any script in `bin`.

---

## 3. Sequence alignment & comparison

Dependency-free local alignment, wrappers around external aligners, and the
score/HSP/windowed-identity machinery used to describe and compare alignments.

### `src/dark/local_align.py`

A from-scratch Smith-Waterman local-alignment implementation for two
sequences (used where you want a dependency-free local alignment
instead of shelling out to an external aligner like
`edlib`/MAFFT/needle — see `dark.aligners` for those).

- **`LocalAlignment(seq1, seq2, match=1, mismatch=-1, gap=-1, gapExtend=-1, gapExtendDecay=0.0)`** — Constructor takes two `dark.reads.Read`-like objects (uses `.sequence` and `.id`) and scoring parameters (`match` must make sense as positive, `mismatch`/`gap` must be negative, `gapExtend` must not be positive — raises `ValueError` otherwise, as does passing either sequence as empty). `gapExtendDecay` lets the per-position gap-extension penalty shrink the longer a gap run gets (a float decay added to `gapExtend` for each extra position in the run, capped so the effective extension penalty never becomes positive).
  - **`createAlignment(resultFormat=dict)`** — The main entry point: runs the DP fill-and-traceback and returns either `None` (no alignment found — the whole table collapsed to zero) or a `dict` with keys `cigar` (a CIGAR-like string using `=`/`X`/`I`/`D`, *not* standard SAM CIGAR — note it uses `=`/`X` always, not `M`), `sequence1Start`/`sequence1End`, `sequence2Start`/`sequence2End` (1-based-looking indices from the internal DP table), and `text` (a 3-line `[seq1Line, matchLine, seq2Line]` representation). If `resultFormat=str`, returns a single formatted multi-line string instead. Internal/private helpers (`_initialise`, `_fillAndTraceback`, `_cigarString`, `_formatAlignment`, `_alignmentToStr`) implement the DP table, traceback, and output formatting.

**Used by:** `local-align.py`.

### `src/dark/aligners.py`

Wraps three different ways of aligning sequences — the external tools
MAFFT and EMBOSS needle/stretcher (via `dark.process.Executor`,
i.e. shelling out), and the `edlib` Python binding for a fast,
dependency-light pairwise alignment — behind a common interface, plus
post-processing helpers to clean up `edlib`'s output. This is the
module that `bin/compare-sequences.py`, `bin/compare-aa-sequences.py`,
`bin/fasta-identity-table.py`, `bin/combine-sequences.py`, and
`bin/fasta-diff.py` build on for their "run an aligner" option.

- **`align(reads, aligner, alignerOptions=None, verbose=False, threads=<cpu_count>)`** — The main dispatch entry point: given `aligner` of `"mafft"`, `"needle"`, or `"edlib"`, calls the corresponding function below with default options substituted when `alignerOptions is None`, and returns the resulting `Reads`.
- **`mafft(reads, verbose=False, options=None, threads=None, executor=None, dryRun=False)`** — Writes `reads` to a temp FASTA file, shells out to `mafft`, and reads the aligned result back as a `Reads` instance. If `dryRun=True`, no command is actually executed and an empty `Reads` is returned.
- **`needle(reads, verbose=False, options=None, executor=None, dryRun=False)`** — Aligns exactly the first two of `reads` using EMBOSS `needle` (global Needleman-Wunsch); if `needle` fails specifically because the sequences are too big, transparently falls back to running EMBOSS `stretcher` instead.
- **`removeFirstUnnecessaryGaps(seq1, seq2, gapSymbol='-')`** — Finds the first stretch in two equal-length, already-aligned strings where one side has a gap-run immediately followed by non-gap characters that could be "slid" out to eliminate the gap run without making the alignment worse. Returns `(True, newSeq1, newSeq2)` or `(False, "", "")`.
- **`removeUnnecessaryGaps(seq1, seq2, gapSymbol='-')`** — Repeatedly applies `removeFirstUnnecessaryGaps` until no more improving gap removals can be found.
- **`canonicalizeAllGaps(seq1, seq2)`** — Post-processes an `edlib`-produced alignment to merge many short gap runs into fewer, longer ones so the result looks more like what MAFFT would produce.
- **`edlibAlign(reads, gapSymbol='-', minimizeGaps=True, onlyTwoSequences=True, matchAmbiguous=True, canonicalizeGaps=True)`** — Aligns the first two of `reads` using `edlib`'s global (`"NW"`) alignment mode, with ambiguous-nucleotide-code equalities taken into account if `matchAmbiguous=True`. Raises `ValueError` if more than two reads are passed while `onlyTwoSequences=True`, or if either input sequence already contains `gapSymbol`.
- Module constants: `MAFFT_DEFAULT_ARGS`, `NEEDLE_DEFAULT_ARGS`, `EDLIB_AMBIGUOUS`.

**Used by:** `fasta-identity-table.py`, `combine-sequences.py`, `fasta-diff.py`, `compare-sequences.py`.

### `src/dark/windowedIdentity.py`

Computes nucleotide identity between a reference sequence and a set of
other (aligned, equal-length) sequences in sliding windows along the
genome — the computational core behind
`bin/plot-windowed-identity.py`.

- **`WindowedIdentity(reads: Iterable[DNARead])`** — Constructor just stores the (already-aligned) reads; validation is deferred to first use.
  - **`getIdentity(referenceRegex, window, minWindow, jump, includeRegex='.', sort=False, strict=True, startOffset=0, stopOffset=None)`** — The main method. Picks exactly one read whose id matches `referenceRegex` as the reference (raises `ValueError` if zero or more than one match), picks reads matching `includeRegex` to compare against it, then slides a window of size `window` (advancing by `jump`, never going below `minWindow` at the end) along the reference. For each window and comparison read, computes identity via `dark.dna.compareDNAReads`. Returns `(reference, windowStarts, identity)` where `identity[read]` is the per-window identity-fraction list. Raises `ValueError` for many invalid parameter combinations; all reads must be the same length (checked internally).
- **`addCommandLineOptions(parser)`** — Adds a standard set of argparse options (`--reference`, `--include`, `--window`, `--minWindow`, `--jump`, `--sort`, `--startOffset`, `--stopOffset`) mirroring `getIdentity`'s parameters.

**Used by:** `plot-windowed-identity.py`.

### `src/dark/hsp.py`

Defines the data-holding classes for a single matched region
("high-scoring pair") between a read and a subject sequence, as
produced by alignment tools (BLAST, DIAMOND) and consumed throughout
the `dark.alignments`/`dark.titles`/`dark.blast`/`dark.diamond` code.

- **`_Base`** — Private base class holding the common HSP/LSP fields: `readStart`, `readEnd`, `readStartInSubject`, `readEndInSubject`, `readFrame`, `subjectStart`, `subjectEnd`, `subjectFrame`, `readMatchedSequence`, `subjectMatchedSequence`, `identicalCount`, `percentIdentical`, `positiveCount`, `percentPositive`. Comparison is delegated to `self.score`.
  - `betterThan(score)` — delegates to the score object's own `betterThan`.
  - `toDict()` — returns a plain `dict` snapshot of the fields above.
- **`HSP(score, **kwargs)`** — `_Base` subclass for "numerically higher is better" scores (bit scores); wraps `score` in `dark.score.HigherIsBetterScore`.
- **`LSP(score, **kwargs)`** — `_Base` subclass for "numerically lower is better" scores (e-values); wraps `score` in `dark.score.LowerIsBetterScore`.

**Used by:** Not used directly by any script in `bin` (used pervasively via `dark.alignments`, `dark.titles`, `dark.blast`, `dark.diamond`).

### `src/dark/score.py`

Tiny wrapper classes that make "is this score better than that score"
comparisons polarity-agnostic, so code elsewhere (e.g. `dark.hsp`,
`dark.titles`) doesn't need to know whether higher-is-better (bit
scores) or lower-is-better (e-values) semantics apply.

- **`HigherIsBetterScore(score: float)`** — `<`/`==` compare numerically as normal; `betterThan(score)` is `True` if this instance's score is numerically greater.
- **`LowerIsBetterScore(score: float)`** — Same idea, inverted. Both classes only compare against instances of their own class (`TypeError` if you mix them).

**Used by:** Not used directly by any script in `bin` (used via `dark.hsp.HSP`/`LSP`).

### `src/dark/distance.py`

A single general-purpose string-distance function.

- **`levenshtein(source, target)`** — Computes the `int` Levenshtein (edit) distance between two strings, using a numpy-vectorized, two-row (not full-matrix) dynamic-programming implementation for memory efficiency.

**Used by:** `adaptor-distances.py`, `split-fasta-by-adaptors.py`.

---

## 4. Alignment-result containers, filtering & summarization

The application-agnostic layer above individual BLAST/DIAMOND HSPs: generic
per-read/per-title match containers, filtering, and the shared
read-filtering/editing argparse helpers.

### `src/dark/filter.py`

General-purpose filtering primitives for sequence titles and
sequences, plus two pairs of `argparse` helper functions (mirroring
the `addFASTACommandLineOptions`/`parseFASTACommandLineOptions`
pattern in `dark.reads`) that let any `bin/` script add a large,
standard set of "select/omit" and "edit" options for FASTA/FASTQ
reads. Used by both the BLAST and DIAMOND alignment-filtering code
(via `dark.alignments`) and directly by several `bin/` scripts that
filter/edit raw reads.

- **`TitleFilter(whitelist=None, blacklist=None, whitelistFile=None, blacklistFile=None, positiveRegex=None, negativeRegex=None, truncateAfter=None)`** — Decides whether a sequence *title* is acceptable.
  - **`accept(title)`** — Returns one of three class constants: `REJECT`, `WHITELIST_ACCEPT`, or `DEFAULT_ACCEPT`. Checks, in order: exact whitelist, exact blacklist, a positive regex, a negative regex, and (if `truncateAfter`) de-duplication by truncated title.
- **`SequenceFilter(whitelist=None, blacklist=None, whitelistFile=None, blacklistFile=None, positiveRegex=None, negativeRegex=None)`** — Same idea as `TitleFilter` but applied to the sequence *string* itself.
  - **`accept(sequence)`** — Same three-way return convention.
- **`ReadSetFilter(minNew)`** — Decides whether a new title's set of matching read ids is sufficiently different from previously-seen titles' read sets (e.g. to drop near-duplicate hits).
  - **`accept(title, titleAlignments)`** — `True` if the fraction of *new* read ids is at least `minNew`; records which titles were invalidated otherwise. Raises `AssertionError` if the same title is passed twice.
  - **`invalidates(title)`** — list of other titles rejected because their read sets overlapped too much with `title`'s.
- **`addFASTAFilteringCommandLineOptions(parser)`** — Adds a large set of read-*selection* options: `--minLength`/`--maxLength`, `--maxNFraction`, title/sequence whitelist/blacklist (inline or file), regex filters, `--keepSequences`/`--removeSequences`, `--head`, deduplication options, `--randomSubset`/`--trueLength`, `--sampleFraction`, `--sequenceNumbersFile`.
- **`parseFASTAFilteringCommandLineOptions(args, reads)`** — Reads those option values and calls `reads.filter(...)`, returning the filtered `Reads`.
- **`addFASTAEditingCommandLineOptions(parser)`** — Adds options that *change* reads: site keep/remove ranges, `--removeGaps`, `--truncateTitlesAfter`, `--removeDescriptions`, `--idLambda`/`--readLambda` (each `eval`'d), `--reverse`/`--reverseComplement`, `--rotate`, `--setSite`, upper/lower-case options, `--setQuality`.
- **`parseFASTAEditingCommandLineOptions(args, reads)`** — The parse half of the editing options, returning the edited `Reads`.

**Used by:** `fasta-identity-table.py`, `filter-fasta.py`, `genome-protein-summary.py`, `filter-sam.py`, `sam-coverage.py`, `sam-to-fasta-alignment.py`.

### `src/dark/alignments.py`

The application-agnostic core of dark matter's alignment-handling
machinery: generic containers for one subject match (`Alignment`) and
all of a read's matches (`ReadAlignments`), a generic
`ReadsAlignments` base class that BLAST- and DIAMOND-specific
subclasses extend, a composable filter (`ReadsAlignmentsFilter`) that
implements most of the actual filtering logic for matches, and two
small gap/offset-mapping helpers used when working with aligned
(gapped) sequences. Nothing here is specific to BLAST or DIAMOND — see
`dark.blast.alignments.BlastReadsAlignments` and
`dark.diamond.alignments.DiamondReadsAlignments` for the concrete
subclasses.

- **`Alignment(subjectLength, subjectTitle)`** — Holds one subject sequence's length/title plus a growing list of HSPs/LSPs matched against it.
  - **`addHsp(hsp)`** — Appends an `HSP` or `LSP`.
- **`ReadAlignments(read, alignments=None)`** — A `list` subclass holding all the `Alignment`s found for one query `read`.
- **`bestAlignment(readAlignments)`** — Returns the `Alignment` whose *first* HSP has the best score.
- **`ReadsAlignmentsParams(application, applicationParams=None, subjectIsNucleotides=True, scoreTitle="Score")`** — A plain data holder recording which application produced a set of alignments and how.
- **`ReadsAlignmentsFilter(limit=None, maxAlignmentsPerRead=None, minSequenceLen=None, maxSequenceLen=None, minStart=None, maxStop=None, oneAlignmentPerRead=False, maxHspsPerHit=None, scoreCutoff=None, percentageIdenticalCutoff=None, percentagePositiveCutoff=None, whitelist=None, blacklist=None, whitelistFile=None, blacklistFile=None, titleRegex=None, negativeTitleRegex=None, truncateTitlesAfter=None, taxonomy=None, readIdRegex=None)`** — Builds a `dark.filter.TitleFilter` and/or `dark.taxonomy.LineageFetcher` internally as needed.
  - **`filter(readAlignments)`** — Applies, in a fixed order, a result-count limit, per-read alignment cap, read-id regex, title filter, sequence-length bounds, taxonomic lineage membership, `oneAlignmentPerRead`, `maxHspsPerHit`, score/percent cutoffs, and start/stop bounds. Returns the narrowed `readAlignments`, or `False` if everything was filtered out.
  - **`close()`** — Closes the internal `LineageFetcher`, if one was created.
- **`ReadsAlignments(reads, params, scoreClass=HigherIsBetterScore)`** — Base class for "a whole alignment-results file (or set of files)."
  - **`getSubjectSequence(title)`** — Abstract; raises `NotImplementedError` here.
  - **`hsps()`** — Generator yielding every HSP/LSP across every (filtered) alignment.
  - **`__iter__()`** — Runs each `readAlignments` through every filter added via `self.filter(...)`.
  - **`iter()`** — Placeholder; subclasses override to actually read alignment data.
  - **`filter(**kwargs)`** — Constructs a `ReadsAlignmentsFilter(**kwargs)` and chains it. Returns `self`.
- **`getGappedOffsets(s, gapChar="-")`** — Maps each offset in the *ungapped* version of `s` to the corresponding offset in the gapped string `s`.
- **`alignmentEnd(s, startOffset, length, gapChar="-")`** — Starting at `startOffset` in gapped sequence `s`, scans forward until `length` non-gap characters have been seen and returns the resulting (gapped) index.

**Used by:** Not used directly by any script in `bin` (scripts go
through the BLAST/DIAMOND-specific subclasses instead).

### `src/dark/titles.py`

Groups and summarizes read alignments (HSPs) by the subject ("title")
they matched — the layer above individual `HSP`s and below the overall
`ReadsAlignments` abstraction, used for producing the per-subject
summaries and filtering seen in the BLAST/DIAMOND alignment-panel
scripts.

- **`titleCounts(readsAlignments)`** — Returns a `dict` mapping subject title → `int` count of how many times that title was matched.
- **`TitleAlignment(read, hsps)`** — Pairs one `Read` with the `list[HSP]` it matched a particular title with.
- **`TitleAlignments(subjectTitle, subjectLength)`** — A `list` subclass holding all the `TitleAlignment`s for one subject title, plus subject metadata.
  - **`addAlignment(alignment)`**, **`reads()`/`readCount()`/`readIds()`**, **`hsps()`/`hspCount()`**, **`bestHsp()`/`worstHsp()`**, **`hasScoreBetterThan(score)`**, **`medianScore()`**.
  - **`coverage()`/`coverageCounts()`** — fraction of the subject covered, or per-position coverage counts, via `dark.intervals.ReadIntervals`.
  - **`coverageInfo()`** — `dict` keyed by subject offset, a pileup-like list of `(score, base)` pairs.
  - **`residueCounts(convertCaseTo='upper')`** — `dict` keyed by subject offset, each value a `Counter` of residue→count.
  - **`summary()`** / **`toDict()`** — compact and full serializations.
- **`TitlesAlignments(readsAlignments, scoreClass=None, readSetFilter=None)`** — A `dict` subclass (title → `TitleAlignments`) built from a `ReadsAlignments`-like object.
  - **`filter(minMatchingReads=None, maxMatchingReads=None, minMedianScore=None, withScoreBetterThan=None, minNewReads=None, minCoverage=None, maxTitles=None, sortOn='maxScore', titleRegex=None, negativeTitleRegex=None)`** — mutates and returns `self` with titles removed that fail any given criterion; `minNewReads` uses a lazily-created `dark.filter.ReadSetFilter`; `maxTitles` keeps only the best-scoring titles.
  - **`sortTitles(by)`** — sorts titles by `'length'`, `'maxScore'`, `'medianScore'`, `'readCount'`, or `'title'`.
  - **`hsps()`**, **`summary(sortOn=None)`**, **`tabSeparatedSummary(sortOn=None)`** (fixed field order — downstream code relies on it), **`toDict()`**.

**Used by:** `alignment-panel-civ.py`, `filter-hits-to-fasta.py`, `get-reads.py`, `noninteractive-alignment-panel.py`.

### `src/dark/intervals.py`

Tracks the (possibly overlapping) set of regions on a subject sequence
that are covered by matching reads, and derives coverage statistics
and plot-offset adjustments from that. Used by `dark.titles` (for
`coverage()`/`coverageCounts()`) and by graphics code that wants to
compress large unmatched regions when plotting (`OffsetAdjuster`).

- **`ReadIntervals(targetLength: int)`** — Holds raw `(start, end)` intervals added via `add`.
  - **`add(start, end)`**, **`walk()`** (generator yielding `(type, (start, end))` covering the entire span with no gaps; `type` is `FULL` or `EMPTY`), **`coverage()`** (fraction covered), **`coverageCounts()`** (per-base depth `Counter`).
- **`OffsetAdjuster(intervals=None, base=2.0)`** — Precomputes, from a `ReadIntervals`' empty stretches, offset adjustments that compress those stretches logarithmically rather than linearly.
  - **`adjustments()`**, **`adjustOffset(offset)`**, **`adjustHSP(hsp)`** (mutates an HSP/LSP *in place*).

**Used by:** Not used directly by any script in `bin` (used via `dark.titles.TitleAlignments.coverage`/`coverageCounts`, and by `dark.graphics`/`dark.orfs`-style plotting code for `OffsetAdjuster`).

---

## 5. BLAST

Dark matter's own compact JSON encoding of BLAST results, and the
BLAST-specific subclass of the generic alignments machinery.

### `src/dark/blast/records.py`

A tiny debugging helper for dumping a raw BioPython BLAST record to stdout.

- **`printBlastRecord(record)`** — Prints every top-level attribute of a `Bio.Blast.Record.Blast` instance (except `alignments`/`descriptions`/`reference`), then per-alignment description/alignment fields, then each HSP via `dark.blast.hsp.printHSP`.

**Used by:** Not used directly by any script in `bin`.

### `src/dark/blast/hacks.py`

A single one-off, admittedly hacky function kept around from a
specific past debugging session; not part of the supported API.

- **`printBlastRecordForDerek(record)`** — Skips records with no alignments or whose best e-value is worse than `1e-120`; otherwise prints query name, each alignment's e-value, and (for the first HSP only) sense/half-genome info using a hard-coded coordinate. The docstring says it's called by `bin/print-blast-xml-for-derek.py`, but no such script currently exists in `bin`.

**Used by:** Not used directly by any script in `bin` (the script named in its own docstring does not currently exist).

### `src/dark/blast/params.py`

A small consistency check for comparing the global BLAST parameters
recorded in one results file against those from another, used when
concatenating/combining multiple BLAST runs.

- **`VARIABLE_PARAMS`** — parameter names expected to legitimately differ between files (and so excluded from comparison).
- **`checkCompatibleParams(initialParams, laterParams)`** — Returns `None` if compatible, else a human-readable summary of every difference found.

**Used by:** Not used directly by any script in `bin` (used internally by `dark.blast.alignments.BlastReadsAlignments`).

### `src/dark/blast/score.py`

Pure-math conversions between BLAST bit scores and e-values, following
NCBI's "BLAST Scoring Parameters" reference.

- **`bitScoreToEValue(bitScore, dbSize, dbSequenceCount, queryLength, lengthAdjustment)`** — Converts a bit score to an e-value via an "effective" database size.
- **`eValueToBitScore(eValue, dbSize, dbSequenceCount, queryLength, lengthAdjustment)`** — The inverse.

**Used by:** `bit-score-to-e-value.py`, `e-value-to-bit-score.py`.

### `src/dark/blast/hsp.py`

A debug-print helper plus the function that converts a raw (dict-form)
BLAST HSP into the normalized, 0-based, Python-slicing-friendly
offsets used throughout the rest of dark matter.

- **`printHSP(hsp, indent="")`** — Prints a fixed set of HSP dict fields, one per line.
- **`normalizeHSP(hsp, readLen, blastApplication)`** — Takes a raw BLAST HSP (1-based offsets, a `frame` 2-tuple) and returns a `dict` with 0-based `readStart`/`readEnd`/`readStartInSubject`/`readEndInSubject`/`subjectStart`/`subjectEnd`, all in the subject's coordinate system. Handles reverse-strand offset flipping, `blastx` nucleotide→protein coordinate conversion, and read-overhang computation. Raises `AssertionError` (after dumping diagnostic info) on internal consistency failures.

**Used by:** Not used directly by any script in `bin` (used internally by `dark.blast.conversion.JSONRecordsReader`).

### `src/dark/blast/alignments.py`

The BLAST-specific concrete subclass of
`dark.alignments.ReadsAlignments`: reads one or more BLAST result
files (dark matter's converted JSON) and makes their HSPs available as
filterable per-read alignment collections, with support for
e-value-based scoring quirks (zero e-values) and for looking up
subject sequences from several possible sources.

- **`ZERO_EVALUE_UPPER_RANDOM_INCREMENT = 150`** — width of the random range used when spreading out zero-e-value HSPs for plotting.
- **`BlastReadsAlignments(reads, blastFilenames, databaseFilename=None, databaseDirectory=None, sqliteDatabaseFilename=None, scoreClass=HigherIsBetterScore, sortBlastFilenames=True, randomizeZeroEValues=True)`** — `reads` must be the exact `Reads` given to BLAST, in output order. `blastFilenames` must be `.json`/`.json.bz2` (dark matter's format, produced by `bin/convert-blast-xml-to-json.py`) — raw XML is not accepted here. `databaseFilename`/`sqliteDatabaseFilename` are mutually-exclusive subject-lookup sources; if neither given, falls back to `blastdbcmd` via `dark.ncbidb.getSequence`.
  - **`iter()`** — yields one `ReadAlignments` per read.
  - **`getSubjectSequence(title)`** — looks up and caches the subject sequence, returning an `AARead` or `DNARead` depending on the BLAST application.
  - **`adjustHspsForPlotting(titleAlignments)`** / **`adjustPlot(readsAx)`** — e-value → `-log10` conversion and zero-e-value placeholder handling, for plotting.

**Used by:** `alignments-per-read.py`, `alignment-panel-civ.py`, `check-fasta-json-blast-consistency.py`, `filter-hits-to-fasta.py`, `get-reads.py`, `filter-reads-alignments.py`, `noninteractive-alignment-panel.py`.

### `src/dark/blast/conversion.py`

The two-way conversion between raw BLAST XML output and dark matter's
own compact JSON encoding of BLAST results (one JSON object per line:
a parameters line, then one record-dict line per query).

- **`XMLRecordsReader(filename)`** — Reads BLAST XML (`-outfmt 5`) via BioPython's `Bio.Blast.NCBIXML.parse`.
  - **`records()`** — generator yielding raw BioPython records; sets `self.params` from the first record.
  - **`saveAsJSON(fp)`** — writes the JSON form: a params line, then one line per record (keeping only the fields dark matter actually uses).
- **`JSONRecordsReader(filename, scoreClass=HigherIsBetterScore)`** — Reads dark matter's JSON-encoded BLAST output; `.bz2` files are transparently decompressed. Raises `ValueError` if the file is empty, malformed, or in an old format missing `"application"`.
  - **`readAlignments(reads)`** — generator yielding one `ReadAlignments` per line, matched in lockstep against `reads`; calls `dark.blast.hsp.normalizeHSP` for each HSP. Raises `ValueError` on a query/read id mismatch or a premature end of the reads iterator.

**Used by:** `convert-blast-xml-to-json.py` (`XMLRecordsReader`), `read-blast-json.py` (`JSONRecordsReader`).

---

## 6. DIAMOND

The DIAMOND-specific counterparts of the BLAST machinery above: tabular
output parsing, JSON conversion, DIAMOND-to-SAM conversion, and a small
wrapper for actually running the `diamond` executable.

### `src/dark/diamond/run.py`

A convenience wrapper for actually invoking the external `diamond`
command-line tool from Python — building a temporary DIAMOND database
from a set of subject sequences, searching a set of query reads
against it, and parsing the tabular results back into Python
dicts. Mainly useful for tests and small ad hoc scripts; production
pipelines more commonly run DIAMOND externally and feed its output
through `dark.diamond.conversion`/`dark.diamond.sam` instead.

- **`diamondInstalled()`** — `True` if `diamond help` runs successfully.
- **`DiamondExecutor(dryRun=False)`** — manages a temporary working directory holding subject/query FASTA and DIAMOND's tabular output; usable as a context manager.
  - **`addSubject(subject)`**, **`search(reads, fieldNames=None)`** (rebuilds the database if dirty, runs `diamond blastx --outfmt 6`, yields one `dict` per result row), **`cleanup()`**.
- **`SUBJECTS_FILENAME`, `QUERIES_FILENAME`, `OUTPUT_FILENAME`** — fixed temp-file basenames.

**Used by:** Not used directly by any script in `bin` (used by the test suite to drive real DIAMOND runs).

### `src/dark/diamond/alignments.py`

Holds `DiamondReadsAlignments`, the DIAMOND-specific counterpart to
the BLAST `ReadsAlignments` machinery: it ties a set of query reads to
the (JSON-encoded) DIAMOND hits found for them, and supports the
plotting/filtering code built on top of `ReadsAlignments`.

- **`DiamondReadsAlignments(reads, filenames, databaseFilename=None, databaseDirectory=None, sqliteDatabaseFilename=None, scoreClass=HigherIsBetterScore, sortFilenames=False, randomizeZeroEValues=True)`** — a `ReadsAlignments` subclass. `filenames` are one or more `.json`/`.json.bz2` files produced by `bin/convert-diamond-to-json.py`. Exactly one of `databaseFilename`/`sqliteDatabaseFilename` must be given.
  - **`iter()`** — one `ReadAlignments` per query read.
  - **`getSubjectSequence(title)`** — looks up and caches the `AAReadWithX` for a subject title; raises `KeyError` if not found.
  - **`adjustHspsForPlotting(titleAlignments)`/`adjustPlot(readsAx)`** — e-value-only plotting adjustments (no-ops for bit-score scoring).
- `ZERO_EVALUE_UPPER_RANDOM_INCREMENT = 150`.

**Used by:** `alignment-panel-civ.py`, `filter-reads-alignments.py`, `noninteractive-alignment-panel.py`.

### `src/dark/diamond/hsp.py`

Converts a raw DIAMOND HSP (1-based offsets) into a normalized,
0-based form usable for display/plotting against a subject sequence,
with internal sanity checks.

- **`normalizeHSP(hsp, queryLen, diamondTask)`** — Given an HSP `dict` (1-based `query_start`/`query_end`/`sbjct_start`/`sbjct_end`, `btop`, `frame`) and query length, returns a `dict` with 0-based `readStart`/`readEnd`/`readStartInSubject`/`readEndInSubject`/`subjectStart`/`subjectEnd`. Handles DIAMOND's `blastx` nucleotide→protein coordinate conversion and reversed-frame queries; `readStartInSubject`/`readEndInSubject` can fall outside `[subjectStart, subjectEnd]` if the query overhangs the subject. Raises `AssertionError` (via `_sanityCheck`/`_debugPrint`) on internal inconsistency.

**Used by:** Not used directly by any script in `bin` (used internally by `dark.diamond.conversion.JSONRecordsReader`).

### `src/dark/diamond/sam.py`

Converts DIAMOND tabular (`--outfmt 6`) matches into SAM output — this
is the engine behind `bin/convert-diamond-to-sam.py`. **Note:** that
script was recently flagged (with a startup stderr warning) as
incomplete and not fully trustworthy — treat this module's output with
the same caution until it's been reviewed further.

- **`SAMWriter(ram=False)`** — accumulates SAM body lines (RAM or temp-file-backed) and tracks `@SQ` reference lengths, raising `ValueError` on a conflicting length for the same reference.
  - **`addMatchLine(line, genome, fp)`**, **`save(fp, tmpfp=None, headerLines=None)`**.
- **`_DiamondSAMWriter(genomesProteins, mappingQuality=255, ram=False, keepDescriptions=False)`** — abstract base class; `genomesProteins` is a `dark.civ.proteins.SqliteIndex` instance (built by `bin/make-protein-database.py`).
  - **`_preprocessMatch(diamondStr)`** — parses a tabular line and looks up its genome/protein, raising `ValueError` on a failed lookup.
  - **`_SAMLine(match, protein, genome)`** — builds the TAB-separated SAM line (reverse-complementing as needed, building the CIGAR via `dark.btop.btop2cigar`, converting the protein `sstart` to a nucleotide `POS`, adding an `AS:i:<bitscore>` tag).
  - **`addMatch(diamondStr)`/`save(...)`** — abstract.
- **`SimpleDiamondSAMWriter`** — writes all matches to one combined SAM output.
- **`PerReferenceDiamondSAMWriter(genomesProteins, baseFilenameFunc=None, ..., fpcMaxsize=0)`** — writes one `<basename>.sam` file per matched reference, using a `dark.fpcache.FilePointerCache` to stay under the OS open-file-descriptor limit.

**Used by:** `convert-diamond-to-sam.py`.

### `src/dark/diamond/conversion.py`

Defines dark-matter's DIAMOND tabular ↔ internal-dict ↔ JSON
conversions: reading raw DIAMOND `--outfmt 6` output, converting it to
the project's compact per-line JSON alignment format, and reading that
JSON format back for use as `ReadAlignments`.

- `FIELDS` — the `--outfmt 6` field-name string DIAMOND must be run with for its output to be usable by `bin/convert-diamond-to-sam.py`.
- `DIAMOND_FIELD_CONVERTER` — maps known DIAMOND field names to their value-conversion function.
- **`DiamondTabularFormat(fieldNames=None)`** — generic reader/converter for arbitrary `--outfmt 6` column layouts.
  - **`diamondFieldsToDict(line)`**, **`diamondTabularFormatToDicts(filename)`**.
- **`DiamondTabularFormatReader(filename)`** — reads *raw* DIAMOND `--outfmt 6` output in one of three older fixed column orders (a different convention than `FIELDS` above) and groups consecutive same-query lines into per-query alignment records.
  - **`records()`**, **`saveAsJSON(fp, writeBytes=False)`**.
- **`JSONRecordsReader(filename, scoreClass=HigherIsBetterScore)`** — reads back dark-matter's per-line JSON DIAMOND format, transparently handling `.bz2`.
  - **`readAlignments(reads)`** — one `ReadAlignments` per entry in `reads`, matched by id (tolerating truncated ids).

**Used by:** `convert-diamond-to-json.py`, `convert-diamond-to-sam.py` (for `FIELDS`).

---

## 7. SAM/BAM, read mapping & consensus building

The pysam-based core of SAM/BAM handling, the external-tool wrappers that
produce BAM files in the first place, and the (pure-Python) consensus-calling
engine built on top of it.

### `src/dark/sam.py`

Core SAM/BAM file handling: reference inspection, coverage-depth
computation, a flexible alignment filter (`SAMFilter`), conversion of
CIGAR-aligned reads to padded (gapped) FASTA (`PaddedSAM`), a
reference-to-reference distance matrix built from shared read matches
(`DistanceMatrix`), reference/consensus-building helpers, and summary
statistics over all reads in a SAM/BAM file. Built on top of `pysam`.

- **`SamError`** / **`UnequalReferenceLengthError`** / **`UnknownReference`** / **`UnspecifiedReference`** / **`ReferenceNameMismatchError`** / **`InvalidSAM`** — a small exception hierarchy (all subclass `SamError`).
- **`CONSUMES_QUERY`** / **`CONSUMES_REFERENCE`** — `set`s of `pysam` CIGAR operation codes that consume query/reference bases, used when walking CIGAR strings.
- **`samfile(filename)`** — `@contextmanager` opening a `pysam.AlignmentFile`, guaranteed closed.
- **`samtoolsInstalled()`**, **`samReferences(filenameOrSamfile)`**, **`samReferenceLengths(filenameOrSamfile)`**, **`samReferencesToStr(filenameOrSamfile, indent=0)`**.
- **`coverageDepth(filename, reference=None)`** — per-base coverage depth for one reference (auto-detected if not given, raising `UnspecifiedReference`/`UnknownReference` as appropriate). Returns `(referenceName, depths)`.
- **`SAMFilter`** — the central class for reading and filtering alignments from a SAM/BAM file.
  - `__init__(filename, filterRead=None, referenceIds=None, storeQueryIds=False, dropUnmapped=False, dropSecondary=False, dropSupplementary=False, dropDuplicates=False, keepQCFailures=False, minScore=None, maxScore=None, scoreTag="AS")`.
  - `addFilteringOptions(parser, samfileIsPositional=False, ...)` (`@staticmethod`) — adds the standard `--samfile`/positional, `--referenceId`, drop-flags, score-range options to an `argparse.ArgumentParser`.
  - `parseFilteringOptions(args, filterRead=None, storeQueryIds=False)` (`@classmethod`) — builds a `SAMFilter` from parsed args.
  - `filterAlignment(alignment)` — tests one `pysam` alignment against all configured criteria.
  - `alignments()` — generator of passing alignments; repairs `*`-SEQ secondary/supplementary records per the SAM spec. Raises `InvalidSAM` if no preceding alignment is available to borrow a sequence from.
  - `referenceLengths(getAllReferences=False)`.
- **`PaddedSAM`** — turns SAM/BAM queries into FASTA-style reads padded with gap characters so they all line up against one reference.
  - `__init__(samFilter)` — requires all references to have the same length (raises `UnequalReferenceLengthError` otherwise).
  - `queries(rcSuffix="", rcNeeded=False, padChar="-", queryInsertionChar="N", unknownQualityChar="!", allowDuplicateIds=False, addAlignment=False)` — generator of padded `Read`s, walking the CIGAR to handle matches/insertions/deletions/clipping, handling reverse-strand reads and duplicate-id suffixing.
- **`DistanceMatrix`** — accumulates read/reference match scores from one or more SAM files and computes inter-reference distance/similarity matrices.
  - `addFile(filename, scoreTag=None)`, `score(referenceId, queryId)`, `jaccardDistance(...)`/`soergelDistance(...)`, `matrix(referenceIds=None, metric="soergel", similarity=False, returnDict=False)`, `save(fp)`/`load(fp)`.
- **`getReferenceInfo(bam, bamFilename, bamId=None, referenceFasta=None, fastaId=None, quiet=False)`** — resolves "which reference, and how long is it" for consensus-building, cross-checking a BAM's own reference info against an optional separate reference FASTA. Returns `(bamId, referenceRead_or_None, referenceLength)`.
- **`ReferenceReads`** — per-reference bookkeeping of which reads matched it and how (primary/secondary/supplementary/duplicate/QC-fail, covered offsets); `coverage()`.
- **`ReadsSummary`** — reads a whole SAM/BAM file and builds a `ReferenceReads` per matched reference.
  - `toString(excludeZeroes=False, excludeIfNoAdditional=False, sortBy="count")`, `sortReferences(sortBy="count")`.

**Used by:** `filter-sam.py`, `sam-references.py`, `sam-reference-read-counts.py`, `sam-coverage.py`, `sam-coverage-depth.py`, `plot-coverage-depth.py`, `sam-to-fasta-alignment.py`, `genome-protein-summary.py`, `simple-consensus.py`, `reference-read-scores-to-JSON.py`, `plot-references-by-inter-read-distance.py`.

### `src/dark/consensus.py`

The core (pure-Python) BAM-to-consensus-sequence engine used by
`bin/simple-consensus.py` (via `consensusFromBAM`); handles per-site
base counting, insertions, deletions, and coverage-dependent symbol
choice. **Note:** `bin/simple-consensus.py` was recently flagged (with
a startup stderr warning) as incomplete and not fully trustworthy —
bear that in mind when relying on this module.

- **`ConsensusError`** — raised for consensus-making problems (e.g. an unknown `strategy`).
- **`Insertion(insertionOffset)`** — accumulates per-read insertion data at a single reference offset.
  - `start(anchorOffset)`, `append(base, quality)`, `readCount()`, `resolve()` (returns a `list[dark.dna.Bases]`, one per inserted position).
- **`basesToConsensus(offsetBases, otherBases, originalOffsets, reference, referenceLength, threshold, minCoverage, lowCoverage, noCoverage, deletionSymbol, progress)`** — turns per-site `Bases` tallies into the final consensus `str`.
- **`getConsensusId(bamId, idLambda)`** — builds the consensus id (`f"{bamId}-consensus"` by default, or an `eval`'d user lambda).
- **`consensusFromBAM(bamFilename, bamId=None, referenceFasta=None, fastaId=None, consensusId=None, idLambda=None, threshold=0.8, minCoverage=1, lowCoverage="n", noCoverage="N", deletionSymbol="-", deletionThreshold=0.5, ignoreQuality=False, insertionCountThreshold=5, strategy="fetch", includeSoftClipped=False, compareWithPileupFile=None, progress=False, quiet=False)`** — the main entry point: resolves the reference (via `dark.sam.getReferenceInfo`), collects per-site data, builds and returns a consensus `DNARead`. Raises `UnspecifiedReference`/`UnknownReference`/`UnequalReferenceLengthError`/`ReferenceNameMismatchError` (from `dark.sam`) or `ConsensusError` for an unsupported `strategy` (only `"fetch"` is implemented).
- **`getPairs(bam, bamId, referenceLength, ignoreQuality, includeSoftClipped, progress)`** — builds per-site `correspondences`/`deletions`/`insertions` by delegating to `addPairsInfo` for each read.
- **`addPairsInfo(pairs, cigarOperations, queryId, query, qualities, referenceLength, includeSoftClipped, correspondences, deletions, insertions)`** — per-read helper; raises `ValueError` on mismatched pairs/CIGAR lengths.
- **`fetchConsensus(bam, correspondences, deletions, insertions, reference, referenceLength, noCoverage, deletionThreshold, ignoreQuality, insertionCountThreshold, includeSoftClipped, progress)`** — turns raw per-site data into the offset-indexed structures `basesToConsensus` consumes.
- **`pileupCorrespondences(bam, bamId, referenceLength, includeSoftClipped, progress)`** — an alternative, pysam `pileup`-based way of collecting per-site data, used only for cross-checking.
- **`compareCorrespondences(fp, fetch, pileup, threshold, minCoverage)`** — diagnostic comparison of the two collection strategies, wired up via `compareWithPileupFile`.

**Used by:** `simple-consensus.py`.

### `src/dark/bowtie2.py`

A thin, stateful wrapper for driving the
`bowtie2`/`samtools`/GATK/Picard/bcftools/ivar command-line toolchain
to align, post-process, and call variants from a FASTQ dataset, via
`dark.process.Executor`. All the real work happens in external shell
commands; this class just sequences them and tracks what stage of the
pipeline it's in.

- **`Bowtie2(executor=None, threads=None, verboseFp=None, dryRun=False, reference=None, tempdir=None, tmpChmod=None)`** (`@final`) — creates a private temp directory for intermediate files.
  - `buildIndex(index)` — auto-detects whether `index` is a pre-built bowtie2 index, a FASTA filename, or an NCBI accession, and arranges to use/build accordingly. Must run before `align`.
  - `align(bowtie2Args="--no-unal", fastq1=None, fastq2=None, threads=None, discardSAM=False, readGroup="orig", sampleName="orig")` — raises `ValueError` if `buildIndex` hasn't run or neither FASTQ is given.
  - `makeBAM(...)` / `sort(byName=False)` / `indexBAM()` / `removeDuplicates()` / `markDuplicatesPicard(picardFile)` / `markDuplicatesGATK(threads=None)` / `removePrimers(bedFile)` — successive post-processing steps, several requiring `makeBAM` to have run first.
  - `callHaplotypesGATK(picardJar, vcfFile=None, referenceFasta=None)` / `callHaplotypesBcftools(vcfFile=None, referenceFasta=None)` — variant-calling into a VCF.
  - `isEmpty()`, `mappedAndUnmappedCounts()`, `outputFile()`, `close()`.

**Used by:** `run-bowtie2.py`.

### `src/dark/process.py`

A single utility class for running (or dry-run logging) shell commands
uniformly across the various `bin/run-*`/`bin/make-consensus.py`-style
scripts, with consistent error reporting and an append-only log of
everything that was (or would have been) executed.

- **`Executor(dryRun=False, stdout=None, stderr=None)`** (`@final`) — `self.log` accumulates a human-readable record of every command run.
  - `execute(command, dryRun=None, useStderr=True, stdout=False, stderr=False, **kwargs)` — runs `command` via `subprocess.run(..., check=True)`, logging it either way; in dry-run mode returns `None` without executing. On `CalledProcessError`, optionally prints the error/output to stderr (per `useStderr`) before re-raising.

**Used by:** `compare-aa-sequences.py`, `compare-consensuses.py`, `run-bowtie2.py`, `run-bwa.py`, `extract-reference-from-fasta-and-bam.py`, `make-consensus.py`.

### `src/dark/extract_reference.py`

Implements the two halves of
`bin/extract-reference-from-fasta-and-bam.py`: pulling a single named
reference sequence out of a FASTA file, and filtering a SAM/BAM file
down to just the records aligned against it.

- **`extract_fasta(fasta, ref_id, out_dir, index_fasta, quiet, force, keep_description, e)`** — finds `ref_id` in `fasta`, writes it to `<out_dir>/<ref_id>.fasta`, optionally indexes it with `samtools faidx`, returns its length. Calls `sys.exit(...)` (not an exception) on a pre-existing output file (unless `force`) or an unfindable reference.
- **`extract_sam(bam, ref_id, reference_len, out_dir, keep_sam, force, e)`** — converts to SAM if needed, copies the matching `@SQ`/body lines, sorts/indexes back to BAM. Returns the matching-read count; `sys.exit`s on zero matches, pre-existing outputs, or an unrecognized file suffix.

**Used by:** `extract-reference-from-fasta-and-bam.py`.

### `src/dark/analyze_reads.py`

Small helper module backing `bin/trim-reads.py`: finds a prefix/suffix
common to every read in a FASTA file (e.g. leftover primer/adaptor
sequence) and trims it off.

- **`getPrefixAndSuffix(file_handle)`** — returns `(prefixLength, suffixLength)`: the longest prefix/suffix common to *every* read in the file. Returns `(0, 0)` if any two reads share none, or if there are no reads.
- **`trimReads(prefix, suffix, file_handle)`** — generator yielding each read with `prefix`/`suffix` bases removed.
- `_longestPrefixOfTwoSeqs(a, b)` — internal helper, not public API.

**Used by:** `trim-reads.py`.

### `src/dark/fpcache.py`

A small LFU (least-frequently-used) file-pointer cache, letting code
hold many "simultaneously open" files for writing/reading without
exceeding the OS's per-process open-file-descriptor limit — files that
fall out of the cache are transparently closed, and reopened (in a
different mode) if touched again.

- **`FilePointerCache`** — the public wrapper class.
  - `__init__(maxsize=32, openArgs=None, reopenArgs=None)` — `openArgs`/`reopenArgs` control the kwargs used for the first open vs. a reopen after cache eviction.
  - `open(filename)` — returns an open file pointer, opening/evicting as needed.
  - `close()` — closes every currently-open file. Usable as a context manager.

**Used by:** Not used directly by any script in `bin` (used internally by `dark.diamond.sam.PerReferenceDiamondSAMWriter`, reached via `bin/convert-diamond-to-sam.py`'s `--perReferenceOutput` option).

### `src/dark/progress.py`

A one-function wrapper around the `progressbar` package that makes an
optional, terminal-aware progress bar easy to use with a single `with`
statement.

- **`maybeProgressBar(show, maxValue, prefix)`** — `@contextmanager`. If `show` and stderr is a TTY, yields a configured `progressbar.ProgressBar`; otherwise yields a no-op stand-in with a do-nothing `update` method, so callers never need to branch.

**Used by:** Not used directly by any script in `bin`.

---

## 8. Protein/genome databases, GenBank & taxonomy

Two parallel sets of "protein/genome database" code exist in this codebase —
the plain `dark.proteins`/`dark.genomes` modules (built around NCBI-style
`[pathogen name]`-bracketed titles) and the `dark.civ.proteins` module (built
around dark matter's own "civ"-prefixed sequence id format and its own
SQLite schema). They are **not** interchangeable; several `bin` scripts'
import bugs (see below) have come from confusing the two.

### `src/dark/proteins.py`

Groups DIAMOND/BLAST protein-match results by the pathogen each
matched protein comes from, and renders that grouping as plain text or
as a browsable HTML report (with optional per-pathogen read-count
panel plots). This is the library behind `proteins-to-pathogens.py`'s
pathogen-centric view of `noninteractive-alignment-panel.py`
output. **Note:** this module is unrelated to
`src/dark/civ/proteins.py` below.

- **`splitNames(names)`** — Splits an NCBI-style id of the form `"Protein name [pathogen name]"` into `(proteinName, pathogenName)` using the *last* `[...]` section.
- **`getPathogenProteinCounts(filenames)`** — Returns a `Counter` of pathogen name → protein count across the given FASTA files.
- **`PathogenSampleFiles`** — Caches, for each (pathogen, sample) pair, a de-duplicated FASTA/FASTQ of the reads that matched that pathogen in that sample, plus sample/pathogen index files.
  - `__init__(proteinGrouper, format_="fasta")`, `add(pathogenName, sampleName)`, `lookup(pathogenName, sampleName)` (raises `KeyError` if not yet added), `writeSampleIndex(fp)`/`writePathogenIndex(fp)`, `pathogenIndex(pathogenName)`.
- **`ProteinGrouper`** — accumulates per-sample protein-match summary lines and groups them by pathogen.
  - `__init__(assetDir="out", sampleName=None, sampleNameRegex=None, format_="fasta", proteinFastaFilenames=None, saveReadLengths=False, titleRegex=None, negativeTitleRegex=None, pathogenDataDir="pathogen-data")`.
  - `addFile(filename, fp)` — parses one summary file's lines; raises `ValueError` on a duplicate pathogen/sample/protein combination.
  - `toStr(title="Summary of pathogens", preamble=None)`, `toHTML(pathogenPanelFilename=None, minProteinFraction=0.0, minProteinCount=0, pathogenType="viral", ...)`. Its `minProteinFraction > 0.0 or minProteinCount > 0` filtering branch previously referenced `self.genomeAccessions`/`self._db` (leftover code copied from `civ/proteins.py`'s similarly-named class, neither attribute existing on this non-civ class) and would raise `AttributeError` if exercised — **fixed** (2026-10-24) to use this class's actual data model, `self.pathogenNames`/`self._pathogenProteinCount`, and verified with both a `minProteinCount` and a `minProteinFraction` end-to-end test (sample kept when it passes the threshold, deleted when it doesn't).
  - `pathogenPanel(filename)` — saves a per-pathogen read-count-vs-sample bar-plot panel.

**Used by:** `proteins-to-pathogens.py` imports `ProteinGrouper`. `subset-protein-database.py` previously imported a nonexistent `SqliteIndex` from this module (`from dark.proteins import SqliteIndex` — this module defines no such class, so the import raised `ImportError`); **fixed** (2026-10-24) to `from dark.civ.proteins import SqliteIndex` (the same stale-import bug pattern found and fixed in `bin/convert-diamond-to-sam.py`) and verified end-to-end against a real database built with `SqliteIndexWriter`.

### `src/dark/genomes.py`

Correlates a genome's annotated proteins (from a civ-style
protein/genome database) with read coverage from one or more SAM
files, to compute per-protein coverage statistics.

- **`GenomeProteinInfo`** — Holds, for one genome accession, its proteins (and nucleotide offset ranges) and which read ids cover which genome offsets.
  - `__init__(genomeAccession, proteinGenomeDB, checkTranslations=True)` — `proteinGenomeDB` must duck-type as a `dark.civ.proteins.SqliteIndex`. Verifies (if `checkTranslations`) that each protein's stored translation matches the genome nucleotides. Raises `dark.errors.NoSuchGenomeError` if the accession isn't found.
  - `addSAM(filename, filterAlignment=None)` — reads a SAM/BAM file and records covering read ids per genome offset and per protein.
  - `proteinCoverageInfo(proteinAccession, minReadOffsetCount=None)` — coverage stats for one protein; raises `KeyError`/`ValueError` as appropriate.
  - `readIdsForAllProteins()`.

**Used by:** `genome-protein-summary.py`.

### `src/dark/genbank.py`

Parses the odd little location-range mini-language BioPython produces
when it reads GenBank flat files (e.g. `"join{[126386:126881](-),
[125941:126232](+)}"`), and extracts genome/CDS feature metadata from
BioPython GenBank records.

- **`GenomeRanges`** — Parses a GenBank/BioPython range string into `(start, stop, forward)` 0-based ranges (merging contiguous ones). Raises `ValueError` on malformed input.
  - `circular(genomeLength)` — `True` if the ranges wrap around the genome's origin.
  - `startInGenome(match)` — converts a DIAMOND match's 1-based protein `sstart` to a 0-based genome nucleotide offset, accounting for spliced (multi-range) proteins.
  - `orientations()`, `distinctRangeCount(genomeLength)`.
- **`getSourceInfo(genome, keys=("host", "note", "organism", "mol_type"), logfp=None, indent="")`** — extracts qualifier values from a GenBank record's `source` feature.
- **`getCDSInfo(genome, feature, keys=(...))`** — extracts and validates a CDS feature's protein id, translation, gene/note/product, parsed `GenomeRanges`, circularity, strand. Returns `None` (with a warning) for CDS features that fail validation — treated as "skip", not a hard error.

**Used by:** `genome-protein-summary.py`.

### `src/dark/database.py`

A single helper, `getDatabaseConnection()`, that opens a MySQL
connection (via `mysql.connector`) to a local `ncbi_taxonomy`
database, using the `DBI_USER`/`USER` and `DBI_PASSWORD` environment
variables for credentials. This is the (older, MySQL-based) taxonomy
backend used by `dark.taxonomy.LineageFetcher`; the newer
`dark.taxonomy.Taxonomy` class uses an sqlite3 database instead.

**Not used directly by any script in `bin`** (used only indirectly, via `dark.taxonomy.LineageFetcher`, which itself isn't imported by anything in `bin` either).

### `src/dark/sqlite3.py`

- **`sqliteConnect(filename)`** — A thin wrapper around `sqlite3.connect` that, on failure, re-raises `sqlite3.OperationalError` with the filename included in the message.

**Not used directly by any script in `bin`** (used internally by several other `dark` modules, e.g. `dark.taxonomy`, `dark.civ.proteins`, `dark.fasta`).

### `src/dark/ncbidb.py`

- **`getSequence(title, db="nt")`** — Given a BLAST-hit-style title, shells out to `blastdbcmd` to fetch the corresponding sequence from a local BLAST database. Requires `blastdbcmd` and the named database to exist locally.

**Not used directly by any script in `bin`**.

### `src/dark/entrez.py`

- **`getSequence(title, db="nucleotide")`** — The network equivalent of `dark.ncbidb.getSequence`: fetches the full GenBank record from NCBI over the network via `Bio.Entrez.efetch`, returning `None` on failure (e.g. offline). Sets `Entrez.email` as NCBI requires; subject to NCBI's rate limits.

**Not used directly by any script in `bin`** (used by `dark.features.FeatureList` as the default sequence fetcher).

### `src/dark/civ/proteins.py`

This is the core of the "civ" protein/genome database machinery: it
defines the SQLite schema and read/write access classes for a database
of genome and protein records (built from GenBank or JSON genome
records), plus the `ProteinGrouper`/`PathogenSampleFiles` classes that
turn per-sample alignment-panel protein-match summaries into grouped,
per-pathogen HTML/text/plot reports. The "civ" sequence id format this
module defines and consumes is
`civ|GENBANK|<protein-accession>|GENBANK|<genome-accession>|<product>
[<pathogen-name>]` (field 2 is the protein accession, field 4 is the
genome accession). This whole module is the backbone that
`make-protein-database.py` writes to and that
`alignment-panel-civ.py`, `describe-protein-database.py`,
`convert-diamond-to-sam.py`, `genome-protein-summary.py`, and
`proteins-to-pathogens-civ.py` all read from.

- **`PathogenSampleFiles`** — a cache/writer of per-(pathogen, sample) de-duplicated FASTA/FASTQ files, used while building a `ProteinGrouper` HTML report.
  - `add(genomeAccession, sampleName)` — gathers every read matching any protein of that pathogen in that sample, de-duplicates by id, writes the combined reads to `pathogen-<genomeAccession>-sample-<index>.<format>`, records the de-duplicated count back into the owning `ProteinGrouper`, and returns the filename (cached on subsequent calls).
  - `lookup(genomeAccession, sampleName)` — cache-only lookup; raises `KeyError` if not yet added.
  - `writeSampleIndex(fp)`.
- **`ProteinGrouper`** — reads alignment-panel protein-summary text output and groups it by pathogen (genome accession) and sample, then renders plain text, HTML, or a per-pathogen read-count plot panel. Constructor requires a `SqliteIndex`-like `proteinGenomeDatabase` and a taxonomy database.
  - `addFile(filename, fp)` — ingests one summary file into `self.genomeAccessions` (a nested `genomeAccession -> sampleName -> {"proteins": {...}, "uniqueReadCount": None}` dict). Raises `ValueError` if a protein can't be found in the database (commonly means the summary was generated against a different database version) or on a duplicate combination.
  - `toStr(title=None, preamble=None, pathogenType="viral")` — plain-text summary (no de-duplicated read counts — those need `toHTML`'s setup).
  - `toHTML(...)` — the main interactive HTML report (collapsible sections, taxonomy links, optional read-count panel image, optional taxonomy tree); internally calls `_computeUniqueReadCounts` first.
  - `pathogenPanel(filename)` — the per-pathogen read-count-vs-sample plot panel image, usable standalone or via `toHTML`.
  - Class constants: `VIRALZONE`/`ICTV` link prefixes, `READCOUNT_MARKER`, `READ_AND_HSP_COUNT_STR_SEP`.
- **`_Genome`** / **`_GenomeLocation`** / **`_GenomeFeature`** — small adapters wrapping a plain `dict` (as loaded from JSON) so it exposes the same attributes a BioPython `SeqIO.parse` GenBank record would, letting the rest of the module treat JSON-sourced and GenBank-sourced genomes identically.
- **`SqliteIndexWriter`** — creates/updates the SQLite database of genomes and proteins (tables `genomes`/`proteins`), simultaneously writing protein FASTA (with "civ"-format ids) to a file pointer. Constructor takes a filename (or `":memory:"`) and an optional `fastaFp` (default `sys.stdout`).
  - `addGenBankFile(filename, taxonomyDatabase, dnaOnly=False, rnaOnly=False, allowedTaxonomicRanks=None, minGenomeLength=None, maxGenomeLength=None, excludeExclusiveHosts=None, excludeFungusOnlyViruses=False, excludePlantOnlyViruses=False, databaseName=None, proteinSource="GENBANK", genomeSource="GENBANK", duplicationPolicy="error", logfp=None)` — parses a GenBank flat file and adds all its genomes/proteins; returns `(genomeSequencesInFile, genomesAdded, proteinsAdded)`. This is what `make-protein-database.py` calls for GenBank-sourced data.
  - `addJSONFile(filename, taxonomyDatabase, ...)` — same, for a single genome supplied as JSON.
  - `addGenome(genome, source, taxonomyId, proteinCount, databaseName, duplicationPolicy="error", logfp=None)` — raises `DatabaseDuplicationError` on a duplicate accession with `"error"` policy.
  - `addProteins(genome, source, ...)` — extracts every CDS feature (via `dark.genbank.getCDSInfo`), writes each translation to the FASTA file with a civ-format id, and calls `addProtein` for each.
  - `addProtein(accession, genomeAccession, sequence, offsets, forward, circular, rangeCount, gene=None, note=None, product=None, duplicationPolicy="error", logfp=None)`.
  - `close()` — creates unique indices, commits, closes. Also a context manager.
  - Class constants: `PROTEIN_ACCESSION_FIELD = 2`, `GENOME_ACCESSION_FIELD = 4`, `SEQUENCE_ID_PREFIX = "civ"`, `SEQUENCE_ID_SEPARATOR = "|"`.
- **`SqliteIndex`** — read-only lookup access to a database built by `SqliteIndexWriter`. This is the class most "civ" `bin` scripts need as `SqliteIndex` (notably `convert-diamond-to-sam.py`, which originally had a stale import pointing at an unrelated same-named class). Constructor takes a filename or an open `sqlite3.Connection`, plus an LRU `lookupCacheSize` (default 1024).
  - `genomeAccession(id_)` / `proteinAccession(id_)` — extract the genome/protein accession from a civ-format id; raise `IndexError` if too few fields.
  - `findGenome(id_)` / `findProtein(id_)` — look up a row (as a `dict`) by full id or accession; `None` if not found (memoized).
  - `findProteinsForGenome(id_)` — generator of all protein rows for a genome.
  - `execute(query, *args)`, `proteinCount()`, `genomeCount()`, `close()` (closes the connection only if this instance opened it itself).

**Used by:** `alignment-panel-civ.py`, `convert-diamond-to-sam.py`, `describe-protein-database.py`, `genome-protein-summary.py`, `make-protein-database.py`, `proteins-to-pathogens-civ.py`.

### `src/dark/taxonomy.py`

Provides lineage/host/classification lookups against a local NCBI
taxonomy database (built by the separate
[`ncbi-taxonomy-database`](https://github.com/acorg/ncbi-taxonomy-database)
project), plus heuristics for classifying viruses as plant-only,
fungus-only, RNA, DNA, or retroviruses, and argparse helpers for
scripts that need a taxonomy database path.

- **`LineageFetcher`** — The older, MySQL-backed lineage lookup, keyed by BLAST-hit gi number; `lineage(title)` walks up the taxonomy tree (cached in memory); `close()`. Superseded by `Taxonomy` below.
- **`Taxonomy`** — The modern, sqlite3-backed taxonomy interface; usable as a context manager.
  - `__init__(dbFilenameOrConnection)`.
  - `lineageFromTaxid(taxid, skipFunc=None, stopFunc=None)` — walks `nodes`/`names` up to the root, returning a `tuple` of `LineageElement(taxid, name, rank)`. Raises `ValueError` on a missing taxid.
  - `lineage(id_, skipFunc=None, stopFunc=None)` *(memoized)* — like above, but `id_` may be an `int` taxid, a `str` accession, or a `str` name.
  - `hostsFromTaxid(taxid)` *(memoized)*, `taxIdFromName(name)`, `hosts(id_)`.
  - `isFungusOnlyVirus(lineage, title=None)` / `isPlantOnlyVirus(lineage, title=None)` *(both memoized)* — heuristic classification combining family/genus/order name lists, a title regex, and (slower) a `hosts()` lookup.
  - `subsetLineageByRanks(lineage, func)` *(static)*, `close()`.
- **`isRetrovirus(lineage)`**, **`isRNAVirus(lineage)`**, **`isDNAVirus(lineage)`** — free functions testing a lineage against fixed name/rank markers.
- **`isAllowedTaxonomicRank(allowedTaxonomicRanks, lineage)`**, **`formatLineage(lineage, namesOnly=False, separator=", ", prefix="")`**, **`lineageTaxonomyLinks(lineage)`**.
- **`Hierarchy`** / **`_HierarchyNode`** — build a lineage tree (rooted at `"Viruses"` by default) and serialize to `jstree`-compatible JSON via `toJSON()`.
- **`addTaxonomyDatabaseCommandLineOptions(parser)`** / **`parseTaxonomyDatabaseCommandLineOptions(args, parser)`** — argparse helper pair adding/resolving `--taxonomyDatabase` (falling back to the `DARK_MATTER_TAXONOMY_DATABASE` environment variable).

**Used by:** `make-protein-database.py`, `get-taxonomy.py`, `proteins-to-pathogens-civ.py`, `get-hosts.py`, `filter-fasta-by-taxonomy.py`.

---

## 9. Phylogenetic trees

### `src/dark/trees.py`

Two independent phylogenetic-tree utilities: one for cleaning up
`iqtree2`-produced Newick trees that mix bootstrap-support and aLRT
values on node labels (using `dendropy`), and one for re-rooting an
`ete4` tree on the branch leading to the MRCA of a given set of tips.

- **`removeAlrt(treeFile, supportThreshold, alrtThreshold, adjustSupportValues, replacementSupportValue, expectedUnlabeledNodeCount, verbose=False, verboseFp=sys.stderr)`** — reads a Newick tree whose internal node labels are `"<alrt>/<support>"` (as `iqtree2` produces when both values are requested — a format only `dendropy` can parse among the tools tried), strips the aLRT part, and (if `adjustSupportValues`) replaces the support label with `replacementSupportValue` on branches whose support meets `supportThreshold` but aLRT falls below `alrtThreshold`. Returns the resulting Newick `str`. Raises `ValueError` on a malformed label or an unexpected count of unlabeled internal nodes (pass `-1` to skip that check).
- **`ete4root(tree, tipNames, detach=False, scale=1.0)`** — re-roots an `ete4` `Tree` on the branch leading into the MRCA of `tipNames` (`sys.exit`s with a clear message if a name matches zero or more than one node); handles an already-trifurcating root by re-rooting on the excluded branch instead, and no-ops (with a stderr note) if the MRCA is already the root. `scale` divides the new root's two child-branch lengths by the given factor. If `detach=True`, the outgroup subtree is detached after rerooting and the remaining subtree is returned instead. `ete3root` is a backward-compatible alias.

**Used by:** `curate-tree-ete.py` (`ete4root`), `remove-alrt-from-tree.py` (`removeAlrt`).

---

## 10. Visualization: plotting & HTML

Two near-duplicate pairs of plotting/HTML code exist here: the general-purpose
`dark.graphics`/`dark.html` and the "civ"-database-specific
`dark.civ.graphics`/`dark.civ.html`. Both build BLAST/DIAMOND alignment-panel
figures and HTML index pages; the civ versions look up subject accession
numbers through a local protein/genome database rather than online.

### `src/dark/civ/graphics.py`

Plotting code for the "civ" pathogen-grouping pipeline: draws
BLAST/DIAMOND alignment graphs and panels of them against a
protein/genome SQLite database (`dark.civ.proteins.SqliteIndex`),
using matplotlib, and can emit an HTML summary via
`dark.civ.html.AlignmentPanelHTMLWriter`. It is a near-duplicate of
`dark.graphics` (see below) specialized to accept a `proteinGenomeDB`
and look up accession numbers through it rather than through online
feature lookup.

- **`report(msg)`** — print a timestamped progress message to stderr.
- **`alignmentGraph(titlesAlignments, title, accession, addQueryLines=True, showFeatures=True, logLinearXAxis=False, logBase=1.1, rankScores=False, createFigure=True, showFigure=True, readsAx=None, imageFile=None, quiet=False, idList=False, xRange="subject", maxBluePlotReads=None)`** — draws a single alignment graph for one matched subject `title`: whisker lines for each query's full extent plus a solid HSP-match line, at score height, optionally colored per-read via `idList`. `logLinearXAxis=True` logarithmically compresses large gap regions on the X axis (via `dark.intervals.OffsetAdjuster`); `maxBluePlotReads` caps drawn lines. Returns a `dict` with `adjustOffset`/`features`/`minX`/`maxX`/`minY`/`maxY` — needed by callers (e.g. `alignmentPanel`) for consistent alignment/scaling across multiple graphs.
- **`alignmentPanel(titlesAlignments, proteinGenomeDB, sortOn="maxScore", idList=False, equalizeXAxes=False, xRange="subject", logLinearXAxis=False, rankScores=False, showFeatures=True, logBase=1.1, subjectType="protein", maxBluePlotReads=None)`** — a 5-column grid of small `alignmentGraph` subplots, one per matched subject title, sorted by `sortOn`. Accession numbers come from `proteinGenomeDB.proteinAccession`/`.genomeAccession` depending on `subjectType`. Shows the figure. Raises `ValueError` for an invalid `xRange`.
- **`alignmentPanelHTML(titlesAlignments, proteinGenomeDB, sortOn="maxScore", outputDir=None, idList=False, equalizeXAxes=False, xRange="subject", logLinearXAxis=False, logBase=1.1, rankScores=False, showFeatures=True, subjectType="protein", maxBluePlotReads=None)`** — like `alignmentPanel`, but saves one PNG per subject into `outputDir` and writes an `index.html` via `AlignmentPanelHTMLWriter`. `outputDir` is required; `equalizeXAxes=True` raises `NotImplementedError`.
- Module constants: `SMALLEST_LOGGED_GAP_TO_DISPLAY = 20`, `Y_AXIS_UPPER_PADDING = 1.01`, `DEFAULT_LOG_LINEAR_X_AXIS_BASE = 1.1`.

**Used by:** `alignment-panel-civ.py`.

### `src/dark/graphics.py`

The general-purpose (non-civ) counterpart of `dark.civ.graphics`:
BLAST/DIAMOND alignment-graph plotting plus a few standalone
sequence-visualization helpers (dot-plot alignment, amino-acid
property/cluster plots). Subject feature/accession information comes
from online lookups (via `dark.features`) rather than a local civ
database, and it additionally supports drawing open reading frames and
per-base-colored query rendering, which `dark.civ.graphics` does not.

- **`report(msg)`** — same timestamped stderr progress printer.
- **`alignmentGraph(titlesAlignments, title, addQueryLines=True, showFeatures=True, logLinearXAxis=False, logBase=1.1, rankScores=False, colorQueryBases=False, createFigure=True, showFigure=True, readsAx=None, imageFile=None, quiet=False, idList=False, xRange="subject", showOrfs=True)`** — same core behavior as the civ version, plus: `showOrfs=True` adds ORF tracks (via `dark.orfs`, auto-disabled for protein subjects); `colorQueryBases=True` renders each query as a per-nucleotide-colored `dark.baseimage.BaseImage` grid instead of a line. No `proteinGenomeDB`/accession parameter — the title itself is used for display and feature lookup.
- **`alignmentPanel(titlesAlignments, sortOn="maxScore", idList=False, equalizeXAxes=False, xRange="subject", logLinearXAxis=False, rankScores=False, showFeatures=True, logBase=1.1)`** — same grid structure as the civ version, without accession lookups (ORFs/colored bases disabled in panel mode).
- **`alignmentPanelHTML(titlesAlignments, sortOn="maxScore", outputDir=None, idList=False, equalizeXAxes=False, xRange="subject", logLinearXAxis=False, logBase=1.1, rankScores=False, showFeatures=True, showOrfs=True)`** — HTML-with-per-subject-PNGs variant using `dark.html.AlignmentPanelHTMLWriter`.
- **`scoreGraph(titlesAlignments, find=None, showTitles=False, figureWidth=5, figureHeight=5)`** — grid of sorted per-title score line plots, optionally highlighting titles matching a `find` predicate. The docstring flags this as "probably bit rotted."
- **`scatterAlign(seq1, seq2, window=7)`** — quick dot-plot alignment visualization: scatter-plots every position pair where the same `window`-sized substring occurs in both sequences.
- **`plotAAProperties(sequence, propertyNames, showLines=True, showFigure=True)`** — plots (via `dark.aa.propertiesForSequence`) per-residue property values; always returns the computed `dict` even if not shown. Raises `ValueError` for an unknown property.
- **`plotAAClusters(sequence, propertyNames, showLines=True, showFigure=True)`** — same, for property *cluster numbers* (via `dark.aa.clustersForSequence`).
- Module constants: `QUERY_COLORS`, `DEFAULT_BASE_COLOR`, `TABLEAU20`, `SMALLEST_LOGGED_GAP_TO_DISPLAY`, `Y_AXIS_UPPER_PADDING`, `DEFAULT_LOG_LINEAR_X_AXIS_BASE`.

**Used by:** `noninteractive-alignment-panel.py`.

### `src/dark/html.py`

HTML/IPython-display rendering helpers for summarizing BLAST/DIAMOND
alignment results — both for Jupyter-notebook inline display and for
the standalone-HTML alignment-panel output written by
`dark.graphics.alignmentPanelHTML`.

- **`NCBISequenceLinkURL(title, field=None, delim="|")`** — extracts a field from a title and builds an NCBI nucleotide-database URL for it. Raises `IndexError` if the field doesn't exist.
- **`NCBISequenceLink(title, field=None, delim="|")`** — same, wrapped in an HTML `<a>` tag.
- **`summarizeTitlesByTitle/ByCount/ByLength/ByMaxScore/ByMedianScore(titlesAlignments, limit=None)`** — each sorts a `dark.titles.TitlesAlignments` by the named attribute and returns an `IPython.display.HTML` listing read count, length, best/median score, and an NCBI link per title.
- **`AlignmentPanelHTMLWriter`** — accumulates one PNG-plus-metadata entry per subject (`addImage`) and, on `close()`, writes `index.html` (thumbnail table + detail sections, each with a matching-reads FASTA/FASTQ link and an optional features link) plus `style.css`.
  - `__init__(outputDir, titlesAlignments)`, `addImage(imageBasename, title, graphInfo)`, `close()`.
- **`readCountText(readCountColors, count, linkText=None)`** — produces a CSS-classed HTML `<span>` for a count per a `dark.colors.ColorsForCounts` instance (or plain text if `readCountColors is None`).

**Used by:** Not used directly by any script in `bin` (imported by `dark.graphics`, which `bin/noninteractive-alignment-panel.py` actually uses).

### `src/dark/civ/html.py`

A small module providing the HTML writer that
`dark.graphics`/`dark.civ.graphics` use to produce the multi-image
"alignment panel" index page (a grid of thumbnail graphs, each linking
to a full-size version with match details, a link to the reads that
produced it, and, where available, feature/NCBI-link annotations).

- **`AlignmentPanelHTMLWriter`** — accumulates one entry per matched subject (via `addImage`), then writes out `index.html` and `style.css` (via `close`). Constructor takes the output directory and a `dark.titles.TitlesAlignments` instance.
  - `addImage(imageBasename, accession, title, graphInfo)`.
  - `close()` — writes `index.html` (thumbnail table, detail sections each with a reads-file link that saves the matching reads as gzipped FASTA/FASTQ, and an optional features-file link) and `style.css`.

**Used by:** Not used directly by any script in `bin` — used internally by `src/dark/graphics.py` and `src/dark/civ/graphics.py`, which in turn back `bin/alignment-panel-civ.py` and `bin/noninteractive-alignment-panel.py`.

### `src/dark/colors.py`

A small utility for mapping integer counts to colors/CSS class names
based on user-specified thresholds, used to color-code read counts in
HTML/plot output (e.g. in `dark.html.readCountText`).

- **`ColorsForCounts`** — holds a sorted (highest threshold first) list of `(threshold, color)` pairs built from an iterable of `"value color"` strings; a `(0, defaultColor)` entry is auto-appended if no zero-threshold entry was given. Raises `ValueError` on malformed entries.
  - `thresholdToCssName(threshold)`, `thresholdForCount(count)`, `colorForCount(count)`.

**Used by:** `proteins-to-pathogens-civ.py`.

### `src/dark/dimension.py`

A single generic multi-dimensional-counter (odometer-style) iterator, used to walk row/column coordinates when laying out grid panels of plots.

- **`dimensionalIterator(dimensions, maxItems=-1, start=None)`** — given a list of `n` positive integers (or `'*'` for unbounded), yields `n`-tuples enumerating every combination in odometer order. `maxItems` caps output (`-1` means unlimited — required if any dimension is `'*'`); `start` resumes from a specific coordinate. Raises `ValueError` for a non-positive dimension or mismatched `start` length.

**Used by:** Not used directly by any script in `bin` (used internally
by `dark.graphics`/`dark.civ.graphics` to lay out panel-of-plots
grids).

### `src/dark/baseimage.py`

A tiny helper for building an RGB pixel grid representing per-base
coloring of aligned sequence data, used by
`dark.graphics.alignmentGraph`'s `colorQueryBases` option to render a
query sequence base-by-base instead of as a single-color line.

- **`BaseImage`** — a 2D array of RGB triples (`numpy`, all-ones/white initially) representing an `xRange` × `yRange` grid, each optionally "blown up" into an `xScale` × `yScale` pixel block. `.data` is the raw array (suitable for `matplotlib.pyplot.imshow`).
  - `set(x, y, value)` — color the block for grid position `(x, y)`.

**Used by:** Not used directly by any script in `bin` (used internally by `dark.graphics.alignmentGraph`).

### `src/dark/features.py`

Fetches and plots GenBank sequence features (CDS, mat_peptide, etc.)
alongside BLAST alignment-panel figures, annotating a matplotlib
figure with colored bars for each feature and a legend describing
them.

- **`Feature`** — Wraps a single BioPython GenBank feature with 0-based `start`/`end`, an assignable plot `color`, and `legendLabel()` (a human-readable description).
- **`FeatureList(list)`** — fetches a sequence's GenBank record (via a `sequenceFetcher`, default `dark.entrez.getSequence`) and populates itself with `Feature`s for whichever features have a type in `wantedTypes`, colored along a "coolwarm" colormap. Sets `self.offline = True` (and stays empty) if there's no network access.
- **`_FeatureAdder`** — abstract base coordinating "fetch features, then plot them, unless there are too many or we're offline." `add(fig, title, minX, maxX, offsetAdjuster=None, sequenceFetcher=None)`; caps display at `MAX_FEATURES_TO_DISPLAY = 50`.
- **`ProteinFeatureAdder`** / **`NucleotideFeatureAdder`** — concrete subclasses with different `DATABASE`/`WANTED_TYPES` and layouts (nucleotide features laid out by reading frame).

**Not used directly by any script in `bin`** (used internally by `dark.graphics`, which several scripts/notebooks use for alignment-panel plotting).

---

## 11. General utilities

Grab-bag helpers used throughout the rest of the package.

### `src/dark/utils.py`

A grab-bag of generic helpers used across the codebase: file-name
sorting, transparent compressed-file opening, numeric range-string
parsing, formatted count/percentage strings, interval compression, and
a copy of scikit-learn's `entropy` function.

- **`numericallySortFilenames(names)`** — Sorts filenames by the numeric prefix of their basename (so `"10.json.bz2"` sorts after `"1.json.bz2"` but before `"2.json.bz2"` numerically, not lexically).
- **`median(numbers)`** — Like `statistics.median`, but raises `ValueError("arg is an empty sequence")` on empty input.
- **`asHandle(fileNameOrHandle, mode="rt", encoding="UTF-8")`** *(context manager)* — The central "open this, whatever it is" helper: opens a `str`/`Path` (transparently using `gzip.open`/`bz2.open` for `.gz`/`.bz2`), or yields an already-open file-like object unchanged (so it is *not* closed on exit in that case).
- **`parseRangeString(s, convertToZeroBased=False, asList=False)`** — Parses `"1-5,12,100-200"` into a `set` (or ordered `list`) of `int`s. Raises `ValueError` on malformed input.
- **`parseRangeExpression(s, convertToZeroBased=False)`** — A more powerful sibling also understanding `&`, `|`, and `()` around range groups, by substituting ranges with Python `set` literals and `eval`-ing the result. Returns `set()` for empty input. **Security note:** uses `eval` on input-derived text — safe for its intended shape, not for untrusted input.
- **`parseStartEndLabelColor(regionStr)`** — Parses a 1-based `"start:end[:label[:color]]"` region spec into `(start, end, label, color)` with `start` converted to 0-based. Raises `ValueError` on invalid bounds.
- **`baseCountsToStr(counts)`** / **`nucleotidesToStr(nucleotides, prefix="")`** — Format base→count data as human-readable strings.
- **`countPrint(mesg, count, len1, len2=None)`** — the formatting helper behind `dark.dna.matchToString`'s per-field output lines (`"<mesg>: <count>"` or `"<mesg>: <count>/<len1> (<pct>%)"`, two-length variant if `len1 != len2`).
- **`pct(a, b, commas=False, places=3, denom=True)`** — Formats `"<a>/<b> (<pct>%)"`. Asserts `0 <= a <= b`.
- **`cd(newdir)`** *(context manager)* — Temporarily `chdir`s, restoring on exit.
- **`take(iterator, n)`** — generator of successive (up to) `n`-item chunks.
- **`readLabels(fp)`** — Reads a `name<TAB>replacement` file into a `dict` (for Newick tree relabeling).
- **`matchOffset(leftPaddedReference, leftPaddedQuery)`** — the 0-based offset into the *unpadded* reference where the query's real content begins.
- **`openOr(filename, mode="r", defaultFp=None, specialCaseHyphen=True)`** *(context manager)* — yields `defaultFp` if `filename` is `None` (or `"-"`, if `specialCaseHyphen`), else opens `filename`.
- **`intsToIntervals(numbers)`** — sorted `list` of `(start, end)` tuples for each maximal run of consecutive integers.
- **`intsToStringIntervals(numbers, sep="-", minimize=False)`** — `str` form (`"5"`/`"5-9"`); `minimize=True` further abbreviates same-length ranges (e.g. `"2016-20"`). This is what produces the compressed gap/no-coverage location lists in `dark.dna.matchToString`'s output.
- **`gmt()`** — current UTC time as `"YYYY-MM-DD HH:MM:SS (UTC)"`.
- **`entropy(labels)`** — a vendored copy of scikit-learn's now-renamed/deprecated entropy helper (natural log). **`entropy2(labels)`** — the same, base-2.
- **`assignRegionRows(regions)`** — greedy interval-scheduling layout: assigns each `(start, end, label, color)` region to the lowest-numbered non-overlapping row, for stacking plot annotations.

**Used by:** `alignment-panel-civ.py`, `compare-sequences.py`, `compare-aa-sequences.py`, `fasta-split.py`, `msa-find-and-extract.py`, `genome-protein-summary.py`, `plot-coverage-depth.py`, `plot-references-by-inter-read-distance.py`, `sam-coverage-depth.py`, `sam-match-lengths.py`, `noninteractive-alignment-panel.py`, `sam-to-fasta-alignment.py`.

### `src/dark/errors.py`

Four tiny `Exception` subclasses used as specific error types elsewhere in the package: `DatabaseDuplicationError`, `NoSuchGenomeError`, `NoSuchProteinError`, `ReadLengthsNotIdenticalError`. Each is a one-line docstring with no additional behavior.

**Used by:** `fasta-variable-sites.py`, `genome-protein-summary.py`.

### `src/dark/__init__.py`

The package root module. Defines `__version__` (read from installed
package metadata, falling back to `"unknown"`), raises an `Exception`
at import time if running on Python older than 3.10, and defines the
`File` type alias (`TextIO | BinaryIO | str | bytes | Path`) used
throughout the codebase's type hints to mean "anything that can be
treated as a readable/writable file: an open handle, a filename, or
raw bytes."

**Used by:** `dark-matter-version.py` (for `__version__`), `fasta-subtraction.py`, `simple-consensus.py` (also for `__version__`).

### `src/dark/ipynb.py`

Not a real module of its own functionality — just a convenience
re-export list, meant to be used as `from dark.ipynb import *` inside
an IPython/Jupyter notebook to pull in the most commonly-needed names
(`FastaReads`, `BlastReadsAlignments`, `TitlesAlignments`,
`titleCounts`, the `summarizeTitlesBy*` family, and the
`alignmentGraph`/`alignmentPanel`/`scoreGraph` plotting functions) in
one line. Defines `__all__` explicitly so the star-import only brings
in those names.

**Not used directly by any script in `bin`** (it's specifically for
interactive notebook use, not the `bin` scripts).

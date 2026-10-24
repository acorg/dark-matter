# The `Read` class hierarchy

This describes the class hierarchy defined in `src/dark/reads.py` that
underlies almost everything else in dark-matter.

Reading FASTA and FASTQ files (to produce instances of the classes
described below) is backed by
[prseq](https://pypi.org/project/prseq/)
([GitHub](https://github.com/virologyCharite/prseq)), a Rust-backed reader
that's roughly twice as fast as an equivalent pure-Python implementation.

## A `Read` is just a record, not necessarily a sequencing read

Despite the name, a `Read` is not required to be an actual sequencing read.
It is really just a record with three attributes:

* `id` — a string identifying the record.
* `sequence` — a string of sequence characters (nucleotides or amino acids).
* `quality` — an optional string of per-character quality scores, which, if
  given, must be the same length as `sequence`.

That's it. The name `Read` is a historical hangover from the early days of
this package, when it was mostly used to process actual sequencing reads
(e.g., from a sequencer, or in a BAM file). These days a `Read` is just as
likely to hold a reference genome, a consensus sequence, a protein, a
primer, or any other named sequence you care to work with — there's no
requirement that it ever came from a sequencing machine. If you're writing
new code and find yourself wanting "a thing with an id and a sequence,"
`Read` (or one of its subclasses, below) is almost certainly what you want,
regardless of whether what you're holding is a "read" in the biological
sense.

## Overview of the hierarchy

```
Read
├── _NucleotideRead          (private; shared DNA/RNA machinery)
│   ├── DNARead
│   │   └── DNAKozakRead
│   └── RNARead
└── AARead
    ├── AAReadWithX
    ├── AAReadORF
    ├── SSAARead
    │   └── SSAAReadWithX
    └── TranslatedRead
```

All of these live in `src/dark/reads.py`. See `doc/library.md` for the full
per-method API reference; this document is about what each class is *for*
and how they relate to one another.

## `Read` — the base class

Holds `id`, `sequence`, and the optional `quality`. Raises `ValueError` at
construction time if a `quality` string is given whose length doesn't match
`sequence`.

On top of just holding those three attributes, `Read` provides:

* Equality, ordering, and hashing based on `(id, sequence, quality)`, so
  reads can be compared, sorted, and put in sets/dicts.
* `len(read)` — the length of `sequence`.
* `read[start:stop]` — slicing, which returns a new instance of the same
  (sub)class, with `quality` sliced the same way if present.
* `toString(format_="fasta")` / `toDict()` / `fromDict(d)` — conversion to
  and from FASTA/FASTQ text or a plain dict (e.g., for JSON serialization).
* `checkAlphabet(count=10)` — checks that the first `count` characters of
  `sequence` are all members of the class's `ALPHABET` (see below), raising
  `ValueError` if not.
* `setSites(...)`, `rotate(...)`, `reverse()`, `find(...)`,
  `getPrefixAndSuffixOffsets(...)`, `walkHSP(...)`, `lowComplexityFraction()`
  — various sequence-editing and inspection methods; see `doc/library.md`
  for details.

`Read.ALPHABET` is an empty `set` by default, meaning "no restriction on
what characters may appear." Subclasses that represent a specific kind of
molecule (DNA, RNA, protein) override it with the set of characters that
are actually valid for that molecule, so that `checkAlphabet` has something
real to check against.

## Nucleotide reads

`_NucleotideRead` is a private base class (not meant to be instantiated
directly) that holds behaviour shared by DNA and RNA: `reverseComplement()`
(returns a new instance of the same class, with `sequence` complemented and
reversed and `quality` reversed) and `translations()` (a generator yielding
all six possible reading-frame translations, as `TranslatedRead` instances
— three forward frames and three reverse-complemented frames).

* **`DNARead`** — `ALPHABET = {A, T, C, G}`. Adds `findORF(offset, ...)`,
  which scans codon-by-codon from a given offset looking for an open
  reading frame, returning a dict describing what it found (length, whether
  a start/stop codon was present, the nucleotide sequence, and its
  translation).
* **`RNARead`** — `ALPHABET = {A, T, C, G, U}`. Otherwise just like
  `DNARead`, but with an RNA-appropriate complement table (so `U` pairs
  with `A` rather than `T` pairing with `A`).
* **`DNAKozakRead`** — a `DNARead` subclass representing one Kozak
  consensus sequence found within some other sequence. Besides the usual
  `id`/`sequence`, it records the `originalRead` it was found in, its
  `start`/`stop` offsets within that read, and a `kozakQuality` score (the
  percentage of the five variable Kozak positions that match the preferred
  consensus base). Instances are produced by `dark.dna.findKozakConsensus`
  — you wouldn't normally construct one by hand.

## Amino acid reads

**`AARead`** is the base for all protein-sequence reads. `ALPHABET` is the
set of standard 1-letter amino acid codes. Its `checkAlphabet` does one
extra thing beyond the base class: since the DNA alphabet (`A`, `C`, `G`,
`T`) is a subset of the amino acid alphabet, it specifically raises
`ValueError` if a long-enough sequence looks like it's actually DNA that
was mistakenly passed to `AARead` — a common and otherwise-silent mistake.
`AARead` also provides `properties()`/`propertyDetails()` (per-residue
physicochemical property lookups, via `dark.aaVars`) and `ORFs(openORFs=False)`,
which scans an already-translated amino acid sequence for `M`...`*`-delimited
open reading frames, yielding `AAReadORF` instances.

* **`AAReadWithX`** — identical to `AARead`, except `ALPHABET` also
  includes `'X'` (the standard code for "unknown amino acid"). Use this
  instead of plain `AARead` whenever your protein sequences might contain
  `X`s (e.g., translations of genomic regions with ambiguous or unknown
  nucleotides).
* **`AAReadORF`** — represents one open reading frame found within a
  larger `AARead` by `AARead.ORFs()`. Besides the ORF's own `id`/`sequence`,
  it records the `start`/`stop` offsets within the original read and
  whether the ORF is `openLeft`/`openRight` — i.e., whether the ORF
  potentially extends past the edge of the original sequence because no
  start/stop codon was found there (as opposed to one genuinely being
  absent). You wouldn't normally construct one of these directly.
* **`SSAARead`** ("secondary-structure AA read") — pairs an amino acid
  `sequence` with a `structure` string of the same length (PDB secondary
  structure annotation, e.g., helix/sheet/coil codes). It has no `quality`
  string. Its `toString()` emits two FASTA records (sequence, then
  structure) in the "fasta-ss" format used by `dark.fasta_ss.SSFastaReads`
  and `--fasta-ss` on the command line.
* **`SSAAReadWithX`** — `SSAARead` with the `'X'`-extended alphabet, the
  same relationship `AAReadWithX` has to `AARead`.
* **`TranslatedRead`** — represents one DNA→AA translation of some
  `originalRead` (a nucleotide read), tagged with which `frame` (0, 1, or
  2) and whether it came from the `reverseComplemented` strand. These are
  what `_NucleotideRead.translations()` yields. `maximumORFLength(openORFs=True)`
  returns the length of the longest ORF found within it (via the inherited
  `AARead.ORFs()`).

## Picking a class by name

Several `bin` scripts let you specify, via a `--readClass` command-line
option, which concrete class should be used to hold the sequences they
read. The mapping from such a name (a plain string, e.g. `"DNARead"`) to
the actual class is `dark.reads.readClassNameToClass`, a `dict` covering
`Read`, `DNARead`, `RNARead`, `AARead`, `AAReadWithX`, `AAReadORF`,
`SSAARead`, `SSAAReadWithX`, and `TranslatedRead`. A related function,
`getUnambiguousBases()`, returns a parallel mapping from the same class
names to the set of *unambiguous* bases/residues for that class (so, for
example, `DNARead` maps to `{A, C, G, T}` without any of the IUPAC
ambiguity codes `dark.dna.AMBIGUOUS` also allows in practice).

## `Read` vs `Reads`

Don't confuse the singular `Read` classes described here with the plural
`Reads` class (and its subclasses `ReadsInRAM`, `FastaReads`, `FastqReads`,
`SSFastaReads`, etc.) — `Reads` is a *collection* of `Read` instances (with
filtering, iteration, and file-reading behaviour of its own), not another
kind of individual record. See `doc/library.md`'s entry for
`src/dark/reads.py` for the full `Reads` API.

# dark-matter

(Pure) Python scripts and libraries for bioinformatics.

Source code: [github.com/acorg/dark-matter](https://github.com/acorg/dark-matter)

## Overview

### Utility scripts

The `bin` directory contains over 100 utility scripts. See
[Scripts](scripts.md) for the full list, grouped into categories such as:

1. FASTA/FASTQ: general sequence inspection, filtering, and manipulation
1. DNA/amino-acid translation and codon tools
1. Pairwise sequence comparison
1. Multiple-sequence-alignment (MSA) and alignment-position analysis
1. Consensus-sequence building
1. Read mapping and SAM/BAM analysis
1. BLAST/DIAMOND hit processing and alignment panels
1. Protein/genome SQLite databases and GenBank data
1. Taxonomy
1. Phylogenetic trees (Newick)

### Library functions

The `dark` library, under `src/dark`, is documented in [Library](library.md)
for developers who want to write their own code against it, grouped into
categories such as:

1. Core reads & sequence I/O
1. DNA & amino-acid biology
1. Sequence alignment & comparison
1. Alignment-result containers, filtering & summarization
1. BLAST
1. DIAMOND
1. SAM/BAM, read mapping & consensus building
1. Protein/genome databases, GenBank & taxonomy
1. Phylogenetic trees
1. Visualization: plotting & HTML
1. General utilities

See [Reads](reads.md) for a dedicated explanation of the central `Read`
class hierarchy that underlies most of the above.

## Installation

```sh
pip install dark-matter
```

or, using [uv](https://docs.astral.sh/uv/):

```sh
uv add dark-matter
```

dark-matter supports Python 3.10 through 3.15.

## License

MIT. See
[LICENSE](https://github.com/acorg/dark-matter/blob/master/LICENSE).

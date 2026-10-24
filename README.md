[![PyPI](https://img.shields.io/pypi/v/dark-matter.svg)](https://pypi.org/project/dark-matter/)
[![Python Version](https://img.shields.io/pypi/pyversions/dark-matter.svg)](https://pypi.org/project/dark-matter/)
[![Downloads](https://img.shields.io/pypi/dm/dark-matter.svg)](https://pypi.org/project/dark-matter/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/acorg/dark-matter/actions/workflows/nox.yml/badge.svg)](https://github.com/acorg/dark-matter/actions/workflows/nox.yml)
[![Docs](https://img.shields.io/badge/docs-acorg.github.io-blue.svg)](https://acorg.github.io/dark-matter/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Last Commit](https://img.shields.io/github/last-commit/acorg/dark-matter)](https://github.com/acorg/dark-matter/commits/master)
[![Open Issues](https://img.shields.io/github/issues/acorg/dark-matter)](https://github.com/acorg/dark-matter/issues)

# Dark matter

(Pure) Python scripts and libraries for bioinformatics.

## Overview

### Utility scripts

The `bin` directory contains over 100 utility scripts. These fall into
the following categories:

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

## Installation

```sh
pip install dark-matter
```

or, using [uv](https://docs.astral.sh/uv/):

```sh
uv add dark-matter
```

dark-matter supports Python 3.10 through 3.15.

## Documentation

The full documentation is at
[acorg.github.io/dark-matter](https://acorg.github.io/dark-matter/),
built from the Markdown files below:

* [`doc/scripts.md`](doc/scripts.md) — the command-line scripts in `bin`,
  for people who just want to run existing tools.
* [`doc/library.md`](doc/library.md) — the `dark` Python library under
  `src/dark`, for developers who want to write their own code against it.
* [`doc/reads.md`](doc/reads.md) — the central `Read` class hierarchy.

## Developers

```sh
git clone https://github.com/acorg/dark-matter.git
cd dark-matter
uv sync
```

Run the tests with

```sh
make pytest
```

or directly via

```sh
uv run pytest
```

## License

MIT. See [LICENSE](LICENSE).

[![PyPI](https://img.shields.io/pypi/v/dark-matter.svg)](https://pypi.org/project/dark-matter/)
[![Python Version](https://img.shields.io/pypi/pyversions/dark-matter.svg)](https://pypi.org/project/dark-matter/)
[![Downloads](https://img.shields.io/pypi/dm/dark-matter.svg)](https://pypi.org/project/dark-matter/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://app.travis-ci.com/acorg/dark-matter.svg?branch=master)](https://app.travis-ci.com/acorg/dark-matter)

## Dark matter

A collection of Python tools for filtering and visualizing
[Next Generation Sequencing](https://en.wikipedia.org/wiki/DNA_sequencing#Next-generation_methods)
reads.

## Installation

On Linux (at least) you will need to first:

```sh
sudo apt install zlib1g-dev libbz2-dev liblzma-dev
```

Then you should be able to install via

```sh
$ pip install dark-matter
```

Some additional information is available for
[Linux](doc/linux.md), [OS X](doc/mac.md), and [Windows](doc/windows.md).

## Optional dependencies

Not all of the following are mandatory. It depends which part of the
dark-matter code you try to use.

### Picard

`run-bowtie2.py` can use [Picard](https://broadinstitute.github.io/picard/)
to mark duplicates.

### GATK

`run-bowtie2.py` uses [gatk](https://gatk.broadinstitute.org) if you call
it with `--markDuplicatesGATK` or `--callHaplotypesGATK`. If you want to do
this you'll need to download the `gatk` zip file, unzip it, and either put
the directory where you unzip it into your shell's path or else move the
contents of the zip file into a directory already in your path.

### bcftools

You'll need
[bcftools](http://samtools.github.io/bcftools/howtos/index.html) if you
want to make consensuses using the `make-consensus.py` script. Either
follow
[the installation instructions](http://samtools.github.io/bcftools/howtos/install.html)
or, if you use [brew](https://brew.sh/) (or Linux brew), `brew install bcftools`.

### samtools

The `run-bowtie2.py` and `run-bwa.py` scripts both make use of
[samtools](https://www.htslib.org/).

### EMBOSS

The `compare-sequences.py` script (if called with `--align`) requires the
`needle` and (possibly) `stretcher` from
[EMBOSS](http://emboss.sourceforge.net/).

### iPython Notebook

If you are using dark matter in an
[iPython Notebook](https://ipython.org/notebook.html), you should install
the requirements in `requirements-ipynb.txt`.

## Development

If you are using dark matter as a developer, you should install the
requirements in `requirements-dev.txt`.

#!/usr/bin/env python

from setuptools import setup


# Modified from http://stackoverflow.com/questions/2058802/
# how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package
def version():
    import os
    import re

    init = os.path.join("dark", "__init__.py")
    with open(init) as fp:
        initData = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]+)['\"]", initData, re.M)
    if match:
        return match.group(1)
    else:
        raise RuntimeError("Unable to find version string in %r." % init)


# Explicitly list bin scripts to be installed, seeing as I have a few local
# bin files that are not (yet) part of the distribution.
scripts = [
    "bin/aa-info.py",
    "bin/aa-to-dna.py",
    "bin/aa-to-properties.py",
    "bin/adaptor-distances.py",
    "bin/alignment-panel-civ.py",
    "bin/alignments-per-read.py",
    "bin/bit-score-to-e-value.py",
    "bin/cat-json-blast-records.py",
    "bin/check-fasta-json-blast-consistency.py",
    "bin/codon-distance.py",
    "bin/compare-consensuses.py",
    "bin/compare-sequences.py",
    "bin/convert-blast-xml-to-json.py",
    "bin/convert-diamond-to-json.py",
    "bin/convert-diamond-to-sam.py",
    "bin/convert-sam-to-fastq.sh",
    "bin/create-newick-relabeling-output.py",
    "bin/dark-matter-version.py",
    "bin/describe-protein-database.py",
    "bin/dna-to-aa.py",
    "bin/download-genbank.sh",
    "bin/download-refseq-viral-fasta.sh",
    "bin/download-refseq-viral-gbff.sh",
    "bin/e-value-to-bit-score.py",
    "bin/extract-ORFs.py",
    "bin/fasta-base-indices.py",
    "bin/fasta-count.py",
    "bin/fasta-count-chars.py",
    "bin/fasta-coverage.py",
    "bin/fasta-diff.sh",
    "bin/fasta-find.py",
    "bin/fasta-identity-table.py",
    "bin/fasta-ids.py",
    "bin/fasta-join.py",
    "bin/fasta-lengths.py",
    "bin/fasta-match-offsets.py",
    "bin/fasta-sequences.py",
    "bin/fasta-sort.py",
    "bin/fasta-split-by-length.py",
    "bin/fasta-split-by-id.py",
    "bin/fasta-split.py",
    "bin/fasta-subset.py",
    "bin/fasta-subtraction.py",
    "bin/fasta-to-phylip.py",
    "bin/fasta-translate.py",
    "bin/fasta-variable-sites.py",
    "bin/fastq-set-quality.py",
    "bin/filter-fasta-by-complexity.py",
    "bin/filter-fasta-by-taxonomy.py",
    "bin/filter-fasta.py",
    "bin/filter-hits-to-fasta.py",
    "bin/filter-reads-alignments.py",
    "bin/filter-sam.py",
    "bin/format-fasta.py",
    "bin/genbank-grep.py",
    "bin/genome-protein-summary.py",
    "bin/get-features.py",
    "bin/get-hosts.py",
    "bin/get-reads.py",
    "bin/get-taxonomy.py",
    "bin/graph-evalues.py",
    "bin/ids.py",
    "bin/local-align.py",
    "bin/make-consensus.py",
    "bin/make-fasta-database.py",
    "bin/make-protein-database.py",
    "bin/ncbi-fetch-id.py",
    "bin/newick-to-ascii.py",
    "bin/noninteractive-alignment-panel.py",
    "bin/parse-genbank-flat-file.py",
    "bin/plot-references-by-inter-read-distance.py",
    "bin/position-summary.py",
    "bin/pre-commit.sh",
    "bin/print-read-lengths.py",
    "bin/proteins-to-pathogens.py",
    "bin/proteins-to-pathogens-civ.py",
    "bin/randomize-fasta.py",
    "bin/read-blast-json.py",
    "bin/read-blast-xml.py",
    "bin/reference-read-scores-to-JSON.py",
    "bin/relabel-newick-tree.py",
    "bin/reverse-complement.py",
    "bin/run-bwa.py",
    "bin/run-bowtie2.py",
    "bin/sam-coverage.py",
    "bin/sam-coverage-depth.py",
    "bin/sam-to-fasta-alignment.py",
    "bin/sam-reference-read-counts.py",
    "bin/sam-references.py",
    "bin/sff-to-fastq.py",
    "bin/simple-consensus.py",
    "bin/split-fasta-by-adaptors.py",
    "bin/subset-protein-database.py",
    "bin/summarize-fasta-bases.py",
    "bin/summarize-reads.py",
    "bin/translate.py",
    "bin/trim-primers.py",
    "bin/trim-reads.py",
    "bin/window-split-alignment.py",
    "bin/write-htcondor-job-spec.py",
]

setup(
    name="dark-matter",
    version=version(),
    packages=["dark", "dark.blast", "dark.diamond", "dark.civ"],
    url="https://github.com/acorg/dark-matter",
    download_url="https://github.com/acorg/dark-matter",
    author="Terry Jones, Barbara Mühlemann, Tali Veith, Sophie Mathias, Udo Gieraths, Nikolai Zaki",
    author_email="tcj25@cam.ac.uk",
    keywords=["virus discovery"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    license="MIT",
    description=(
        "Python classes and utility scripts for working with genetic sequence data."
    ),
    scripts=scripts,
    install_requires=[
        "biopython>=1.71",
        "bz2file>=0.98",
        "cachetools>=3.1.0",
        "edlib>=1.3.9",
        "Cython>=0.29.16",
        "ipython>=3.1.0",
        "matplotlib>=1.4.3",
        "mysql-connector-python>=8.1.0",
        "numpy>=1.14.2",
        "progressbar2>=4.0.0",
        "pysam>=0.15.2",
        "pyzmq>=14.3.1",
        "requests>=2.18.4",
        "simplejson>=3.5.3",
        "types-cachetools>=5.5.0",
        "types-requests>=2.32.0",
    ],
)

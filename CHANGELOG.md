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

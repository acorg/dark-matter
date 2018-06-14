## 3.0.15 June 14, 2018
Made a change in `dark/proteins.py`, to make the `minProteinFraction` work on
a per sample basis, not per pathogen.

## 3.0.14 June 12, 2018
Fixed a bug in `diamond/alignments.py` that crept in after version
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

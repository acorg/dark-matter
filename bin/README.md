This directory contains scripts that make use of the dark matter library
functions.

## filter-fasta.py

`filter-fasta.py` has many uses. Run it with `--help` to see options. Two
that may not be apparent are:

* Running it with no arguments will read FASTA and write the same FASTA,
  but sequences will be printed on a single line. This can be very useful
  if you want to perform some simple modification on a FASTA file. For
  example, you could add `/1` to the name of every sequence as follows:

  `filter-fasta.py < filename.fasta | awk 'NR % 2 == 1 {printf "%s/1\n", $0}'`

  This works because the id lines in the FASTA are on the odd-numbered
  lines. If `filter-fasta.py` were not used, this trick would need to be
  done in a different way.

  This becomes much more useful when dealing with FASTQ files, as one
  cannot simply and reliably differentiate between quality strings and
  lines that have sequences or ids. But with

  `filter-fasta.py --fastq < filename.fastq | awk 'NR % 4 == 1 {...}'`

  this becomes trivial.
* It is also easy to convert FASTQ to FASTA:

  `filter-fasta.py --fastq --saveAs fasta < filename.fastq`

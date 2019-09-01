# THIS IS ALL OUT OF DATE! IGNORE.

# Making a protein database from the RVDB nucleotide database

There is a `Makefile` in the `misc` directory that can perform most of the
steps described below for you. To use it, try something like:

```sh
$ mkdir ~/protein-database
$ cd ~/protein-database
$ ln -s ~/dark-matter/misc/Makefile-protein-database Makefile
$ make download
$ make
```

and cross your fingers.  Else, read on and do these things one at a time.

## Executables

There are scripts run below called `download-genbank.sh` and
`make-protein-database.py`. They live in the
[dark-matter](https://github.com/acorg/dark-matter) `bin` directory, which
should be in your `PATH`. Try `pip install dark-matter` if you don't
already have it installed.

## Get the clustered RVDB nucleotide database

Get `C-RVDBv15.1.fasta` (or a later version) from
[here](https://hive.biochemistry.gwu.edu/rvdb).

`C-RVDBv15.1.fasta` contains 699,223 viral nucleotide genomes. We want to
make a database of the proteins in the sequences (actually, the coding
regions).

This is an alternative to using the protein version of the RVDB (available
[here](https://rvdb-prot.pasteur.fr/)) which doesn't do exactly what we
need (we need more intermediate information, such as exactly where the
protein sequences fall in the nucleotide genomes).

## Download the GenBank nucleotide to accession to taxonomy id file

You will need `nucl_gb.accession2taxid.gz` which you can download from
[here](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)
(there's no need to uncompress it).

### Download the GenBank info for all nucleotide genomes

The following uses GNU `parallel` to download files containing 1000
[GenBank flat-file](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)
format records for the accession numbers of the sequences in
`C-RVDBv15.1.fasta`:

```sh
$ download-genbank.sh ids C-RVDBv15.1.fasta accession-to-name
```

If you have an NCBI API key, put its value into your environment in
`NCBI_API_KEY`. If that variable is set the script will send 10 requests
per second to NCBI (instead of 3). See
[here](https://www.ncbi.nlm.nih.gov/books/NBK25497/) for how to obtain an
NCBI API key.


This creates a directory called `ids` (specified by the 3rd argument,
above), holding files with names like

```
$ ls -1 ids/ids_aaaa*
ids/ids_aaaaa.genbank
ids/ids_aaaab.genbank
ids/ids_aaaac.genbank
ids/ids_aaaad.genbank
ids/ids_aaaae.genbank
```

### Accession numbers to names

The above also creates a file called `accession-to-name` (specified by its
3rd argument) that maps accession numbers to nucleotide sequence names (as
in `C-RVDBv15.1.fasta`). That file will be needed below when we make the
protein file for DIAMOND and put the nucleotide sequence name into square
brackets at the end of the line.

## Make the protein FASTA and database

Run

```sh
$ make-protein-database.py --gb ids/ids_*.genbank \
      --databaseFile genomes-proteins.db \
      --nucleotideAccessionToTaxidFile nucl_gb.accession2taxid.gz \
      --accessionToNameFile accession-to-name \
      --force > proteins.fasta
```

This creates `proteins.fasta` which contains an amino acid sequence for all
coding regions in all the nucleotide sequences in `C-RVDBv15.1.fasta`, as
extracted from the `ids/ids_*.genbank` files downloaded above.

The `proteins.fasta` file can be given to `DIAMOND` to make its database.

It also writes a `genomes-proteins.db` file. This is an sqlite3 database with
tables as follows:

```
CREATE TABLE proteins (
    accession VARCHAR UNIQUE PRIMARY KEY,
    sequence VARCHAR NOT NULL,
    offsets VARCHAR NOT NULL,
    reversed INTEGER,
    gene VARCHAR,
    note VARCHAR
);

CREATE TABLE genomes (
    accession VARCHAR UNIQUE PRIMARY KEY,
    length INTEGER,
    proteinCount INTEGER,
    taxid INTEGER
);
```

The former holds information about the protein sequences that will be used
later, principally for determining coverage of the nucleotide sequence the
proteins came from.  The latter contains the number of proteins (coding
regions, actually) for a given nucleotide genome.

## Make a nucleotide sequence id to taxonomy id database

Later, when we match (using DIAMOND) protein sequences, we will be able to
extract their corresponding nucleotide genome accession number from their
sequence id (it's the 2nd accession number in a name like
`acc|GENBANK|YP_009137150.1|GENBANK|NC_001798|neurovirulence protein`).
We'll want to look up the corresponding NCBI taxonomy id, which we can do
either by the protein accession number or the nucleotide one (since they're
obviously from the same organism). Because there are fewer nucleotide
sequences we here make a small database mapping nucleotide accession number
to taxid.

So:

```sh
$ ./add-accession-taxids.py
```

which adds taxonomy ids to all the accession numbers in the `genomes`
database table above. This is used in later processing and filtering.

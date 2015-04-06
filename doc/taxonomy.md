## Installing NCBI taxonomy databases

In order to be able to filter by taxonomic level, you need make a MySQL database
with the NCBI taxonomy information, as below.

### Download the taxonomy database files

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)
[ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz)
These files, one for nucleotide sequences and one for proteins, contain
rows that map the gi number of each NCBI sequence to an NCBI taxonomy id.

If you only care about one of the two, don't download the other, and skip
the corresponding `LOAD` command in the "Load taxonomy data" section below.

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
This `tar` file contains multiple files, the ones we'll need are nodes.dmp and names.dmp.

### Uncompress / extract the needed files

```sh
$ gunzip gi_taxid_nucl.dmp.gz gi_taxid_prot.dmp.gz
$ tar xfvz taxdump.gz
```

### Install MySQL

If you already have MySQL installed, skip this step and the creation of the
`ncbi_taxonomy` database, and go to "Create tables" below.

#### On Mac OS X using brew

```sh
$ brew install mysql
$ unset TMPDIR
$ mysql_install_db --verbose --user=`whoami` --basedir="$(brew --prefix mysql)" --datadir=/usr/local/var/mysql --tmpdir=/tmp
$ mysql.server start
$ /usr/local/opt/mysql/bin/mysql_secure_installation
```

#### On Windows

See: [mySQL installer](http://dev.mysql.com/downloads/windows/installer/).

### Create a database

```sh
$ mysql -uroot -p
mysql> CREATE DATABASE ncbi_taxonomy;
mysql> USE ncbi_taxonomy
```

#### Create tables

```
mysql> DROP TABLE IF EXISTS gi_taxid, names, nodes;
mysql> CREATE TABLE gi_taxid (gi INT, taxID INT);
mysql> CREATE TABLE names (taxID INT, divider1 VARCHAR(300), name VARCHAR(300), divider2 VARCHAR(300), unique_name VARCHAR(300), divider3 VARCHAR(300), name_class VARCHAR(300));
mysql> CREATE TABLE nodes (taxID INT, divider1 VARCHAR(300), parent_taxID INT, divider2 VARCHAR(300), rank VARCHAR(300));
```

#### Load taxonomy data

The following assume you started `mysql` in the directory to the place
where you uncompressed and un-tarred the NCBI taxonomy files. If you
didn't, you'll need to adjust the path names below to refer to the taxonomy
files.

The first two commands below will take some time (possibly hours) to
complete.

If you didn't download both the nucleotide and protein files above, omit
the corresponding `LOAD` command below.

```
mysql> LOAD DATA LOCAL INFILE 'gi_taxid_nucl.dmp' INTO TABLE gi_taxid;
mysql> LOAD DATA LOCAL INFILE 'gi_taxid_prot.dmp' INTO TABLE gi_taxid;
mysql> LOAD DATA LOCAL INFILE 'names.dmp' INTO TABLE names;
mysql> LOAD DATA LOCAL INFILE 'nodes.dmp' INTO TABLE nodes;
```

#### Add indices to the databases

The first command below will take some time (possibly hours) to complete.

```
mysql> ALTER TABLE gi_taxid ADD INDEX (gi);
mysql> ALTER TABLE nodes ADD INDEX (taxID);
mysql> ALTER TABLE names ADD INDEX (taxID);
```

### Clean up

None of the files you downloaded or extracted from the `tar` file are
needed once their contents has been loaded into MySQL.  So unless you have
a reason to keep them, you can remove them and reclaim ~5Gb of disk space.

## Installing NCBI taxonomy databases

In order to be able to filter by taxonomic level, you need make a mysql database
with the NCBI taxonomy information, as below.

### Download the taxonomy database files

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)
This file contains a list that maps the gi number of each database record to a taxonomy id.

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
This directory contains multiple files, the ones we need are nodes.dmp and names.dmp.

### Install mySQL

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
mysql> DROP TABLE IF EXISTS gi_taxid_nucl, names, nodes;
mysql> CREATE TABLE gi_taxid_nucl (gi INT, taxID INT);
mysql> CREATE TABLE names (taxID INT, divider1 VARCHAR(300), name VARCHAR(300), divider2 VARCHAR(300), unique_name VARCHAR(300), divider3 VARCHAR(300), name_class VARCHAR(300));
mysql> CREATE TABLE nodes (taxID INT, divider1 VARCHAR(300), parent_taxID INT, divider2 VARCHAR(300), rank VARCHAR(300));
```

Note that the "DROP TABLE IF EXISTS" command will allow you to upgrade an
existing installation.

#### Load taxonomy data

The following assume you started `mysql` in the directory to the place
where you uncompressed and un-tarred the NCBI taxonomy files. If you
didn't, you'll need to adjust the path names below to refer to the taxonomy
files.

Note that the first command below will take some time to complete.

```
mysql> LOAD DATA LOCAL INFILE 'gi_taxid_nucl.dmp' INTO TABLE gi_taxid_nucl;
mysql> LOAD DATA LOCAL INFILE 'names.dmp' INTO TABLE names;
mysql> LOAD DATA LOCAL INFILE 'nodes.dmp' INTO TABLE nodes;
```

#### Add indices to the databases

Note that the first command below will take some time to complete.

```sh
mysql> ALTER TABLE gi_taxid_nucl ADD INDEX (gi);
mysql> ALTER TABLE nodes ADD INDEX (taxID);
mysql> ALTER TABLE names ADD INDEX (taxID);
```

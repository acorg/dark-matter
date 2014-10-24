## Installing NCBI taxonomy databases

In order to be able to filter by taxonomic level, you need make a mysql database
with the NCBI taxonomy information, as below.

### Download the taxonomy database files

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)
This file contains a list that maps the gi number of each database record to a taxonomy id.

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
This directory contains multiple files, the ones we need are nodes.dmp and names.dmp.

### Install mySQL

#### On Mac OS X using brew

```sh
$ brew install mysql
$ unset TMPDIR
$ mysql_install_db --verbose --user=`whoami` --basedir="$(brew --prefix mysql)" --datadir=/usr/local/var/mysql --tmpdir=/tmp
$ mysql.server start
$ /usr/local/opt/mysql/bin/mysql_secure_installation
```

* make a mySQL database:
```sh
mysql> create database ncbi_taxonomy;
mysql> use ncbi_taxonomy;
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
mysql> CREATE TABLE gi_taxid_nucl (gi INT, taxID INT);
mysql> CREATE TABLE names (taxID INT, divider1 VARCHAR(300), name VARCHAR(300), divider2 VARCHAR(300), unique_name VARCHAR(300), divider3 VARCHAR(300), name_class VARCHAR(300));
mysql> CREATE TABLE nodes (taxID INT, divider1 VARCHAR(300), parent_taxID INT, divider2 VARCHAR(300), rank VARCHAR(300));
```

#### Load taxonomy data

Note that the first command below will take some time to complete.

```
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/gi_taxid_nucl.dmp' INTO TABLE gi_taxid_nucl;
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/gi_taxid_nucl.dmp' INTO TABLE names;
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/nodes.dmp' INTO TABLE nodes;
```

#### Add indices to the databases

Note that the first command below will take some time to complete.

```sh
mysql> ALTER TABLE gi_taxid_nucl ADD INDEX (gi);
mysql> ALTER TABLE nodes ADD INDEX (taxID);
mysql> ALTER TABLE names ADD INDEX (taxID);
```

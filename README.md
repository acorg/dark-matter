## dark-matter

Virus discovery is super fun!

## Installation on Linux

The following works on Ubuntu 14.04 LTS.

### Download & install Freetype2

```sh
$ cd /tmp
$ curl -L 'http://downloads.sourceforge.net/freetype/freetype-2.5.3.tar.bz2' > freetype-2.5.3.tar.bz2
$ tar xfj freetype-2.5.3.tar.bz2
$ cd freetype-2.5.3
$ sed -i  -e "/AUX.*.gxvalid/s@^# @@" -e "/AUX.*.otvalid/s@^# @@" modules.cfg
$ sed -ri -e 's:.*(#.*SUBPIXEL.*) .*:\1:' include/config/ftoption.h
$ ./configure --prefix=/usr --disable-static
$ make
$ sudo make install
```

### Install some required packages

```sh
$ sudo apt-get install python-pip pkg-config python-dev libpng-dev mysql-server libmysqlclient-dev
```

### Install dark matter in a Python virtual environment

```sh
$ sudo pip install virtualenv
$ virtualenv env
$ . env/bin/activate
$ git clone git@github.com:acorg/dark-matter.git
$ cd dark-matter
$ pip install -r requirements.txt
```

### Set PYTHONPATH

You may want to add the `dark-matter` directory (the one that was created
above by `git clone`) to your `PYTHONPATH`.

### Install NCBI taxonomy databases

In order to be able to filter by taxonomic level, you need make a mysql database
with the NCBI taxonomy information according to the following format:

* Download the taxonomy database files from NCBI:

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)
This file contains a list that maps the gi number of each database record to a taxonomy id.

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
This directory contains multiple files, the ones we need are nodes.dmp and names.dmp.

* Install mySQL:
```sh
$ brew install mysql
$ unset TMPDIR
$ mysql_install_db --verbose --user='some username' --basedir="$(brew --prefix mysql)" --datadir=/usr/local/var/mysql --tmpdir=/tmp
$ /usr/local/opt/mysql/bin/mysql_secure_installation
$ mysql.server start
$ /usr/local/opt/mysql/bin/mysql_secure_installation
$ mysql -uroot -p
```

* make a mySQL database:
```sh
mysql> create database ncbi_taxonomy;
```

* create the tables:
```sh
mysql> CREATE TABLE gi_taxid_nucl (gi INT, taxID INT);
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/gi_taxid_nucl.dmp' INTO TABLE gi_taxid_nucl;

mysql> create table names (taxID INT, divider1 VARCHAR(300), name VARCHAR(300), divider2 VARCHAR(300), unique_name VARCHAR(300), divider3 VARCHAR(300), name_class VARCHAR(300));
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/gi_taxid_nucl.dmp' INTO TABLE names;

mysql> create table nodes (taxID INT, divider1 VARCHAR(300), parent_taxID INT, divider2 VARCHAR(300), rank VARCHAR(300));
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/nodes.dmp' INTO TABLE nodes;
```

* Index the databases:
```sh
mysql> ALTER TABLE gi_taxid_nucl ADD INDEX (gi);
mysql> ALTER TABLE nodes ADD INDEX (taxID);
mysql> ALTER TABLE names ADD INDEX (taxID);
```

## Installation on Windows

The following works on Windows 8.1 with Cygwin. Alternative to sudo (which is found in the Linux instructions): run Cygwin as administrator or use `cygstart --action=runas` followed by the command.

### Install pip

See: [pip installation guide](http://pip.readthedocs.org/en/latest/installing.html).

### Install some required Cygwin packages

These can be found and downloaded through the Cygwin setup.exe process or using apt-cyg, if installed:
```
curl
make
python-setuptools
libmysqlclient-devel
mysql
mysqld
libmysqld-devel
libpng-devel
libboost_python-devel
pkg-config
```

### Download & install Freetype2

Ensure you have GNU Make version 3.78.1 or higher by running `make -v`.

```sh
$ cd /tmp
$ curl -L 'http://downloads.sourceforge.net/freetype/freetype-2.5.3.tar.bz2' > freetype-2.5.3.tar.bz2
$ tar xfj freetype-2.5.3.tar.bz2
$ cd freetype-2.5.3
$ sed -i  -e "/AUX.*.gxvalid/s@^# @@" -e "/AUX.*.otvalid/s@^# @@" modules.cfg
$ sed -ri -e 's:.*(#.*SUBPIXEL.*) .*:\1:' include/config/ftoption.h
$ ./configure --prefix=/usr --disable-static
$ make
$ make install
```

### Install dark matter from github

```sh
$ git clone git@github.com:acorg/dark-matter.git
$ cd dark-matter
$ pip install -r requirements.txt
```

### Set PYTHONPATH

You may want to add the `dark-matter` directory (the one that was created
above by `git clone`) to your `PYTHONPATH`.



### Install NCBI taxonomy databases

In order to be able to filter by taxonomic level, you need make a mysql database
with the NCBI taxonomy information according to the following format:

* Download the taxonomy database files from NCBI:

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)
This file contains a list that maps the gi number of each database record to a taxonomy id.

[ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
This directory contains multiple files, the ones we need are nodes.dmp and names.dmp.

* Install mySQL:
See: [mySQL installer](http://dev.mysql.com/downloads/windows/installer/).

* make a mySQL database:
```sh
mysql> CREATE DATABASE ncbi_taxonomy;
mysql> USE ncbi_taxonomy
```

* create the tables:
```sh
mysql> CREATE TABLE gi_taxid_nucl (gi INT, taxID INT);
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/gi_taxid_nucl.dmp' INTO TABLE gi_taxid_nucl;

mysql> create table names (taxID INT, divider1 VARCHAR(300), name VARCHAR(300), divider2 VARCHAR(300), unique_name VARCHAR(300), divider3 VARCHAR(300), name_class VARCHAR(300));
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/gi_taxid_nucl.dmp' INTO TABLE names;

mysql> create table nodes (taxID INT, divider1 VARCHAR(300), parent_taxID INT, divider2 VARCHAR(300), rank VARCHAR(300));
mysql> LOAD DATA LOCAL INFILE 'path-to-dir/nodes.dmp' INTO TABLE nodes;
```

* Index the databases:
```sh
mysql> ALTER TABLE gi_taxid_nucl ADD INDEX (gi);
mysql> ALTER TABLE nodes ADD INDEX (taxID);
mysql> ALTER TABLE names ADD INDEX (taxID);
```

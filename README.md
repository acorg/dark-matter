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


## Installation on Windows

The following works on Windows 8.1 with Cygwin. Alternative to sudo (which is found in the Linux instructions): run Cygwin as administrator or use `cygstart --action=runas` followed by the command.

### Install pip

See: [pip installation guide](http://pip.readthedocs.org/en/latest/installing.html).

### Install some required Cygwin packages

These can be found and downloaded through the Cygwin setup.exe process:
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

### Install dark matter in a Python virtual environment

```sh
$ git clone git@github.com:acorg/dark-matter.git
$ cd dark-matter
$ pip install -r requirements.txt
```

### Set PYTHONPATH

You may want to add the `dark-matter` directory (the one that was created
above by `git clone`) to your `PYTHONPATH`.
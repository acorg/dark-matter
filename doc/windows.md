## Installation on Windows

The following works on Windows 8.1 with [Cygwin](https://www.cygwin.com/).
Alternative to sudo (which is found in the Linux instructions): run Cygwin
as administrator or use `cygstart --action=runas` followed by the command.

### Install pip

See: [pip installation guide](http://pip.readthedocs.org/en/latest/installing.html).

### Install some required Cygwin packages

These can be found and downloaded through the Cygwin `setup.exe` process or
using `apt-cyg`, if installed:

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

### Install dark matter from PyPI

Note that if you are running Python 2, you should use the `requirements-2.txt`
file in the following.

```sh
$ pip install dark-matter
$ cd dark-matter
$ pip install -r requirements-3.txt
```

### Set PYTHONPATH

You may want to add the `dark-matter` directory (the one that was created
above by `git clone`) to your `PYTHONPATH`.

### Install a taxonomy database (optional)

See [taxonomy.md](taxonomy.md) for details.

## Installation on Linux

The following works on Ubuntu 14.04 LTS.

### Install some required packages

```sh
$ sudo apt-get install pkg-config python-dev libpng-dev mysql-server libmysqlclient-dev pkg-config gcc
```

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

### Install dark matter in a Python virtual environment

The dark matter code runs under Python 2 and 3. You need to create a
virtual environment for the Python of your choice and then install the
correct requirements:

#### Python 2

```sh
$ sudo pip install virtualenv
$ virtualenv --python=python2 env
$ . env/bin/activate
$ pip install dark-matter
$ cd dark-matter
$ pip install -r requirements-2.txt
```

#### Python 3

```sh
$ sudo pip install virtualenv
$ virtualenv --python=python3 env
$ . env/bin/activate
$ pip install dark-matter
$ cd dark-matter
$ pip install -r requirements-3.txt
```

### Set PYTHONPATH

You may want to add the `dark-matter` directory (the one that was created
above by `git clone`) to your `PYTHONPATH`.

### Install a taxonomy database (optional)

See [taxonomy.md](taxonomy.md) for details.

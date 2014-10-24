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
$ pip install dark-matter
$ cd dark-matter
$ pip install -r requirements.txt
```

### Set PYTHONPATH

You may want to add the `dark-matter` directory (the one that was created
above by `git clone`) to your `PYTHONPATH`.

### Install a taxonomy database (optional)

See [taxonomy.md](taxonomy.html) for details.

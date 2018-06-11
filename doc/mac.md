## Installation on OS X

The following should works on OS X 10.10.5 (Yosemite) through to 10.13.4
(High Sierra).

# Using Brew and virtualenv

Firstly, install [brew](http://brew.sh/) if you don't have it already.
Then use brew to install
[virtualenv](https://pypi.python.org/pypi/virtualenv).

```sh
$ brew install pyenv-virtualenv
```

## Create a Python virtual environment

```sh
$ virtualenv env
$ . env/bin/activate
```

## Install

Either download a stable version from PyPI using pip:

```sh
$ pip install dark-matter
```

Or clone the dark matter github repo to have the very latest code and
install it manually via `setup.py`.

```sh
$ git clone https://github.com/acorg/dark-matter
$ cd dark-matter
$ python setup.py install
```

# Using pypy

The dark matter code isn't really supported under pypy. Most things should
be fine, but there are currently (2018-06-11) a couple of packages we use
that aren't yet available under pypy. The status is a bit unclear, sorry!

# Install a taxonomy database (optional)

See [taxonomy.md](taxonomy.md) for details.

# Running the tests

If you run the tests using `make check` you may encounter the following
error:

``` 
RuntimeError: Python is not installed as a framework. The Mac OS X backend 
will not be able to function correctly if Python is not installed as a 
framework. See the Python documentation for more information on installing 
Python as a framework on Mac OS X. Please either reinstall Python as a 
framework, or try one of the other backends. If you are using (Ana)Conda 
please install python.app and replace the use of 'python' with 'pythonw'. See
 'Working with Matplotlib on OSX' in the Matplotlib FAQ for more information.
```

You can solve this by editing `~/.matplotlib/matplotlibrc` (you may have to
create the `~/.matplotlib` directory) and inserting the following line:

```
backend: TkAgg
```

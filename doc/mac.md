## Installation on OS X

The following works on OS X 10.10.5 (Yosemite)

# Using Brew and virtualenv

First of all, install [brew](http://brew.sh/) if you don't have it already.
Then use brew to install [virtualenv](https://pypi.python.org/pypi/virtualenv).

```sh
$ brew install pyenv-virtualenv
```

## Create a Python virtual environment

```sh
$ virtualenv env
$ . env/bin/activate
```

## Get the dark matter sources

You can either download a stable (and possibly slightly old) version of
dark matter from PyPI using pip:

Note that if you are using Python 2, use `requirements2.txt` in the
following.

```sh
$ pip install dark-matter
$ cd dark-matter
$ pip install -r requirements.txt
```

Or, clone the dark matter github repo to have the very latest code:

```sh
$ git clone https://github.com/acorg/dark-matter
$ cd dark-matter
$ python setup.py install
```

In this latter case, you may want to add the `dark-matter` directory (the
one created above by `git clone`) to your shell's `PYTHONPATH`.

# Using pypy

Things are a little more awkward to set up using pypy, due to the building
of matplotlib. Here are instructions for brew and virtualenv (first install
those two and make and activate a virtualenv):

```sh
$ git clone https://github.com/acorg/dark-matter
$ cd dark-matter
$ python setup.py install
$ pip install -r requirements-pypy.txt
```

Then install the pypy-specific version of matplotlib:

```sh
$ git clone git@github.com:mattip/matplotlib.git  # This might take a while
$ cd matplotlib
$ python setup.py install
```

# Install a taxonomy database (optional)

See [taxonomy.md](taxonomy.md) for details.

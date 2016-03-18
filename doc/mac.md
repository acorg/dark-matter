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

Note that if you are using Python 2, use `requirements-2.txt` in the
following.

```sh
$ pip install dark-matter
$ cd dark-matter
$ pip install -r requirements-3.txt
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

If you are using pypy 5.0.0 (or later, presumably):

```sh
$ pip install -r requirements-pypy.txt
```

If you're still on pypy 4, comment out the
`git+https://bitbucket.org/pypy/numpy.git` line in `requirements-pypy.txt`
and follow the instructions in that file to install `numpy`. Then run the
`pip install` command above.

# Install a taxonomy database (optional)

See [taxonomy.md](taxonomy.md) for details.

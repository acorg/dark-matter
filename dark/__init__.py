import sys

if sys.version_info < (2, 7):
    raise Exception('The dark matter code needs Python 2.7 or later.')

# Note that the version string must have the following format, otherwise it
# will not be found by the version() function in ../setup.py
__version__ = '1.1.31'

import sys

if sys.version_info < (3, 6):
    raise Exception("The dark matter code needs Python 3.6 or later.")

# Note that the version string below must have the following format,
# otherwise it will not be found by the version() function in ../setup.py
#
# Remember to update ../CHANGELOG.md describing what's new in each version.
__version__ = "4.0.62"

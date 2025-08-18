import sys
from importlib.metadata import PackageNotFoundError, version

if sys.version_info < (3, 10):
    raise Exception("The dark matter code needs Python 3.10 or later.")


try:
    __version__ = version("dark-matter")
except PackageNotFoundError:
    # Package is not installed.
    __version__ = "unknown"

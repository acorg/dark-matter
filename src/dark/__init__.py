import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import BinaryIO, TextIO

if sys.version_info < (3, 10):
    raise Exception("The dark matter code needs Python 3.10 or later.")  # pyright: ignore[reportUnreachable]


try:
    __version__ = version("dark-matter")
except PackageNotFoundError:
    # Package is not installed.
    __version__ = "unknown"


File = TextIO | BinaryIO | str | bytes | Path

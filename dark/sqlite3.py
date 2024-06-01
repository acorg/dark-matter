import sqlite3
from pathlib import Path


def sqliteConnect(filename: str | Path) -> sqlite3.Connection:
    """
    Connect to sqlite3, with a more informative error in case of failure.

    @param : A C{str} or C{Path} database filename.
    @return: An sqlite3 connection.
    """
    try:
        return sqlite3.connect(filename)
    except sqlite3.OperationalError as e:
        # The current (2024-06-01) sqlite3.OperationalError does not
        # include the filename that sqlite3 was trying to open, which
        # is not very helpful. Raise an error that does.
        raise sqlite3.OperationalError(
            f"Could not open sqlite3 database file {str(filename)!r}: {e}"
        )

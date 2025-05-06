from os.path import exists
from threading import RLock
from cachetools import LFUCache, cached


class _FilePointerCache(LFUCache):
    """
    A least-frequently-used file pointer cache class with an eviction method
    that closes open files.
    """

    def popitem(self):
        """
        Evict a cache item.

        @return: The C{str} name of the evicted file.
        """
        _, (filename, fp) = super().popitem()
        fp.close()
        return filename

    def close(self):
        """
        Close all open files in the cache.
        """
        while self.currsize:
            self.popitem()


class FilePointerCache:
    """
    A file pointer cache class for simultaneously opening multiple files,
    without exceeeding the per-process operating system limit on open file
    handles.

    @param maxsize: The C{int} maximum size of the cache.
    @param openArgs: If not C{None}, a C{dict} of keyword arguments to pass to
        open when opening a new file. If C{None}, the file will be opened
        with the default mode ('rt').
    @param reopenArgs: If not C{None}, a C{dict} of keyword arguments to pass
        to open when opening an already existing file. If C{None}, the file
        will be opened with the default mode ('rt').
    """

    def __init__(self, maxsize=32, openArgs=None, reopenArgs=None):
        self._openArgs = openArgs or {}
        self._reopenArgs = reopenArgs or {}
        self._filePointerCache = _FilePointerCache(maxsize=maxsize)
        self._lock = RLock()

        @cached(self._filePointerCache, lock=self._lock)
        def _open(filename):
            if exists(filename):
                fp = open(filename, **self._reopenArgs)
            else:
                fp = open(filename, **self._openArgs)
            return filename, fp

        self._open = _open

    def open(self, filename):
        return self._open(filename)[1]

    def close(self):
        with self._lock:
            self._filePointerCache.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

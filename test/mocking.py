# The following helper function is taken from
# http://www.voidspace.org.uk/python/mock/examples.html#mocking-open

from cStringIO import StringIO
from mock import inPy3k, MagicMock

if inPy3k:
    fileSpec = [
        '_CHUNK_SIZE', '__enter__', '__eq__', '__exit__',
        '__format__', '__ge__', '__gt__', '__hash__', '__iter__', '__le__',
        '__lt__', '__ne__', '__next__', '__repr__', '__str__',
        '_checkClosed', '_checkReadable', '_checkSeekable',
        '_checkWritable', 'buffer', 'close', 'closed', 'detach',
        'encoding', 'errors', 'fileno', 'flush', 'isatty',
        'line_buffering', 'mode', 'name',
        'newlines', 'peek', 'raw', 'read', 'read1', 'readable',
        'readinto', 'readline', 'readlines', 'seek', 'seekable', 'tell',
        'truncate', 'writable', 'write', 'writelines'
    ]
else:
    fileSpec = file


def getFileReader(data):
    content = StringIO(data)

    def reader(n=None):
        if n is None:
            return content.read()
        else:
            return content.read(n)
    return reader


def mockOpen(mock=None, data=''):
    mock = mock or MagicMock(spec=fileSpec)
    handle = MagicMock(spec=fileSpec)
    handle.__enter__.return_value = handle
    handle.write.return_value = None
    handle.read.side_effect = getFileReader(data)
    mock.return_value = handle
    return mock

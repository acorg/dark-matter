from os.path import join, dirname

from cffi import FFI


def fileToString(filename):
    """
    Read a file's contents and return it as a string.

    @param filename: A C{str} file name.
    @return: A C{str} containing the contents of the file.
    """
    with open(join(dirname(__file__), filename)) as f:
        s = f.read()
    return s


ffi = FFI()

header = fileToString('gor4-base.h')

ffi.set_source('dark._gor4', header +
               fileToString('defines.h') +
               fileToString('nrutil.h') +
               fileToString('nrutil.c') +
               fileToString('gor4-base.c') +
               fileToString('api.c'))

ffi.cdef(header)

if __name__ == '__main__':
    ffi.compile()

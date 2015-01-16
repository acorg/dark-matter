cimport cgor4

from os.path import dirname, join
import dark

_DATA_DIR = join(dirname(dark.__file__), '..', 'data', 'gor4')
_SEQUENCES = join(_DATA_DIR, 'New_KS.267.seq')
_SECONDARY = join(_DATA_DIR, 'New_KS.267.obs')

cdef class GOR4:
    cdef cgor4.State *_state
    def __cinit__(self, sequenceFile=_SEQUENCES, secondaryFile=_SECONDARY):
        cdef int error
        self._state = cgor4.initialize(sequenceFile, secondaryFile,
                                       &error)
        if error:
            raise Exception('Error in gor-base.c initialize function.')

        if self._state is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._state is not NULL:
            cgor4.finalize(self._state)

    def __init__(self, sequenceFile=_SEQUENCES, secondaryFile=_SECONDARY):
        pass

    def predict(self, sequence):
        cdef int i
        cdef float *y
        cgor4.predict(self._state, sequence)
        # The gor4-base.c code uses 1-based indexing, unfortunately.
        prob = []
        append = prob.append
        for i in xrange(len(sequence)):
            y = self._state.probai[i + 1]
            append((y[0], y[1], y[2]))
        return {
            'predictions': self._state.predi + 1,
            'probabilities': prob,
        }

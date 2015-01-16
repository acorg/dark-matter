cimport cgor4

cdef class GOR4:
    cdef cgor4.State *_state
    def __cinit__(self, sequenceFilename, secondaryFilename):
        cdef int error
        self._state = cgor4.initialize(sequenceFilename, secondaryFilename,
                                       &error)
        if error:
            raise Exception('Error in gor-base.c initialize function.')

        if self._state is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._state is not NULL:
            cgor4.finalize(self._state)

    def __init__(self, sequenceFilename, secondaryFilename):
        pass

    def predict(self, sequence):
        cdef int i
        cdef float *y
        cgor4.predict(self._state, sequence)
        # The gor-base.c code uses 1-based indexing.
        prob = []
        append = prob.append
        for i in xrange(len(sequence)):
            y = self._state.probai[i + 1]
            append((y[0], y[1], y[2]))
        return self._state.predi + 1, prob

cimport cgor_iv

cdef class GOR_IV:
    cdef cgor_iv.State *_state
    def __cinit__(self, sequenceFilename, secondaryFilename):
        cdef int error
        self._state = cgor_iv.initialize(sequenceFilename, secondaryFilename,
                                         &error)
        if error:
            raise Exception('Error in gor.c initialize function.')

        if self._state is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._state is not NULL:
            cgor_iv.finalize(self._state)

    def __init__(self, sequenceFilename, secondaryFilename):
        pass

    def predict(self, sequence):
        cdef int i
        cdef float *y
        cgor_iv.predict(self._state, sequence)
        # The gor.c code uses 1-based indexing.
        prob = []
        append = prob.append
        for i in xrange(len(sequence)):
            y = self._state.probai[i + 1]
            append((y[0], y[1], y[2]))
        return self._state.predi + 1, prob

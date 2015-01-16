cimport cgor4

from os.path import dirname, join
import dark

_DATA_DIR = join(dirname(dark.__file__), '..', 'data', 'gor4')
_SEQUENCES = join(_DATA_DIR, 'New_KS.267.seq')
_SECONDARY = join(_DATA_DIR, 'New_KS.267.obs')


cdef class GOR4:
    """
    Interface with the C GOR IV code (see src/gor4).

    @param sequenceFile: The C{str} name of a file containing reference
        sequences.
    @param secondaryFile: The C{str} name of a file containing known
        secondary structure information for the sequences in C{sequenceFile}.
    """
    cdef cgor4.State *_state

    def __cinit__(self, sequenceFile=_SEQUENCES, secondaryFile=_SECONDARY):
        cdef int error
        self._state = cgor4.initialize(sequenceFile, secondaryFile,
                                       &error)
        if error:
            raise Exception('Error in gor4-base.c initialize function.')

        if self._state is NULL:
            raise MemoryError()

    def __dealloc__(self):
        """
        Free the storage associated with self._state.
        """
        if self._state is not NULL:
            cgor4.finalize(self._state)

    def __init__(self, sequenceFile=_SEQUENCES, secondaryFile=_SECONDARY):
        pass

    def predict(self, sequence):
        """
        Perform GOR IV prediction on an AA sequence.

        @param sequence: A C{str} sequence of amino acids.
        @return: A C{dict} with 'predictions' and 'probabilities' keys.
            The 'predictions' value is a C{str} of letters from {'H', 'E',
            'C'} for Helix, Beta Strand, Coil.  The probabilities value is
            a C{list} of C{float} triples, one for each amino acid in
            C{sequence}. The C{float} values are the probabilities assigned,
            in order, to Helix, Beta Strand, Coil.
        """
        cdef int i
        cdef float *y
        # The gor4-base.c code uses 1-based indexing, unfortunately. So we need
        # to pad the sequence we pass, and to adjust all received results.
        sequence = 'X' + sequence
        cgor4.predict(self._state, sequence)
        predictions = self._state.predi[1:len(sequence)]
        probabilities = []
        append = probabilities.append
        for i in xrange(1, len(sequence)):
            p = self._state.probai[i]
            append((p[1], p[2], p[3]))
        return {
            'predictions': predictions,
            'probabilities': probabilities,
        }

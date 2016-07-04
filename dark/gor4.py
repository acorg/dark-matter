from ._gor4 import ffi, lib

from os.path import dirname, join
import dark

_DATA_DIR = join(dirname(dark.__file__), 'data', 'gor4')
_SEQUENCES = join(_DATA_DIR, 'new-gor4.seq')
_SECONDARY = join(_DATA_DIR, 'new-gor4.obs')


class GOR4(object):
    """
    Interface with the C GOR IV code (see src/gor4).

    @param sequenceFile: The C{str} name of a file containing reference
        sequences.
    @param secondaryFile: The C{str} name of a file containing known
        secondary structure information for the sequences in C{sequenceFile}.
    """

    def __init__(self, sequenceFile=_SEQUENCES, secondaryFile=_SECONDARY):
        sequenceFile = sequenceFile.encode('UTF-8')
        secondaryFile = secondaryFile.encode('UTF-8')

        error = ffi.new('int *')
        state = lib.initialize(ffi.new('char[]', sequenceFile),
                               ffi.new('char[]', secondaryFile), error)
        if state == 0:
            raise Exception('Error in gor4 initialization.')

        if state == ffi.NULL:
            raise MemoryError()

        # ffi.gc returns a copy of the cdata object which will have the
        # destructor (in this case ``finalize``) called when the Python
        # object is GC'd:
        # https://cffi.readthedocs.org/en/latest/using.html#ffi-interface
        self._state = ffi.gc(state, lib.finalize)

    def predict(self, sequence):
        """
        Perform GOR IV prediction on an AA sequence.

        @param sequence: A C{str} sequence of amino acids.
        @return: A C{dict} with 'predictions' and 'probabilities' keys.
            The 'predictions' value is a C{str} of letters from {'H', 'E',
            'C'} for Helix, Beta Strand, Coil.  The 'probabilities' value is
            a C{list} of C{float} triples, one for each amino acid in
            C{sequence}. The C{float} values are the probabilities assigned,
            in order, to Helix, Beta Strand, Coil.
        """
        # The gor4-base.c code uses 1-based indexing, unfortunately. So we need
        # to pad the sequence we pass, and to adjust all received results.
        sequence = ('X' + sequence).encode('UTF-8')
        lib.predict(self._state, ffi.new('char[]', sequence))
        predictions = ffi.string(self._state.predi[1:len(sequence)])
        if isinstance(predictions, bytes):
            predictions = predictions.decode()
        probabilities = []
        append = probabilities.append
        for i in range(1, len(sequence)):
            p = self._state.probai[i]
            append((p[1], p[2], p[3]))

        return {
            'predictions': predictions,
            'probabilities': probabilities,
        }

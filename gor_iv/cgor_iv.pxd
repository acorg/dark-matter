cdef extern from "gor.h":
    ctypedef struct State:
        char *predi
        float **probai

    State *initialize(char *sequenceFilename, char *secondaryFilename, int *error)
    void finalize(State *state)
    void predict(State *state, char *sequence)

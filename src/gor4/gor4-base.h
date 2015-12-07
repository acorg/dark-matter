typedef struct {
    char **seq;
    char **obs;
    char **title_obs;
    char **title_seq;
    int *temp;
    int *nres;
    int resultLength;
    char *predi;
    float **probai;
} State;

extern State *initialize(char *sequenceFilename, char *secondaryFilename, int *error);
extern void finalize(State *state);
extern void predict(State *state, char *sequence);

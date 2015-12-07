#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _CFFI_
#include "defines.h"
#include "nrutil.h"
#include "gor4-base.h"
#endif

/* Functions in gor.c */
extern void Parameters(int nprot_dbase, int *nres, char **obs, char **seq);
extern void predic(int nres, char *seq, char *pred, float **proba);
extern void First_Pass(int nres, float **proba, char *pred);
extern void Second_Pass(int nres, float **proba, char *pred);

int nprot_dbase = 0;

static int
countProteins(char *filename)
{
    /* Determine the number of proteins in the Kabsch-Sander data base */

    int n = 0;
    FILE *fp;
    char buffer[BUFSIZE];

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Could not open file '%s'.\n", filename);
        return -1;
    }
  
    while (fgets(buffer, BUFSIZE, fp) != NULL) {
        if (buffer[0] == '>' || buffer[0] == '!') {
            n++;
        }
    }

    fclose(fp);
    return n;
}

int
read_file(char *fname, int nprot, char **obs, char **title, int *pnter)
{
    FILE *fp;
    int ip, nres, i;
    int c;
    char *keep;

    fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open file '%s'.\n", fname);
        return 0;
    }

    keep = (char *) malloc((size_t) MAXRES*sizeof(char));

    for(ip = 1; ip <= nprot; ip++) {
        fgets(title[ip], MAXLINE, fp);
        nres = 0;
        while((c = getc(fp)) != '@') {
            if(c == '\n' || c == ' ' || c =='\t'){
                continue;
            }
            nres++;
            if(nres > MAXRES) {
                fprintf(stderr, "The value of MAXRES should be increased: %d", MAXRES);
                return 0;
            }
            if((c >= 'A' && c < 'Z') && c != 'B' && c != 'J' && c != 'O' && c != 'U') {
                keep[nres] = c;
            }
            else {
                fprintf(stderr, "protein: %d residue: %d\n", ip, nres);
                fprintf(stderr, "Invalid amino acid type or secondary structure state: ==>%c<==\n", c);
                return 0;
            }
        }
        while((c = getc(fp)) != '\n'){
            ;
        }
        for(i = 1; i <= nres; i++){
            obs[ip][i] = keep[i];
        }
        pnter[ip] = nres;
    }

    free(keep);
    fclose(fp);
    return 1;
}

static State *
allocate()
{
    State *new = (State *)malloc(sizeof(State));

    new->seq = cmatrix(1, nprot_dbase, 1, MAXRES);
    new->obs = cmatrix(1, nprot_dbase, 1, MAXRES);
    new->title_obs = cmatrix(1, nprot_dbase, 1, MAXLINE);
    new->title_seq = cmatrix(1, nprot_dbase, 1, MAXLINE);
    new->temp = ivector(1, nprot_dbase);
    new->nres = ivector(1, nprot_dbase);
    new->resultLength = 0;
    new->predi = NULL;
    new->probai = NULL;

    return new;
}

State *
initialize(char *sequenceFilename, char *secondaryFilename, int *error)
{
    State *state;

    *error = 0;
    nprot_dbase = countProteins(sequenceFilename);

    if (nprot_dbase == -1){
      *error = 1;
      return 0;
    }

    state = allocate();
    if (!read_file(sequenceFilename, nprot_dbase, state->seq, state->title_seq, state->temp) ||
        !read_file(secondaryFilename, nprot_dbase, state->obs, state->title_obs, state->nres)){
        *error = 1;
        return 0;
    }

    Parameters(nprot_dbase, state->nres, state->obs, state->seq);

    return state;
}

void
predict(State *state, char *sequence)
{
    /*
     * We get passed a sequence that starts at offset one. So we have to
     * subtract one to get its true length.
     */
    int length = (int)strlen(sequence) - 1;

    if (length > state->resultLength){
        /* Allocate more space for the results. */
        if (state->resultLength){
            free(state->predi);
            free(state->probai);
        }
        state->predi = cvector(1, length);
        state->probai = matrix(1, length, 1, 3);
    }

    predic(length, sequence, state->predi, state->probai);
    First_Pass(length, state->probai, state->predi);
    Second_Pass(length, state->probai, state->predi); 
}

void
finalize(State *state)
{
    free_cmatrix(state->seq, 1, nprot_dbase, 1, MAXRES);
    free_cmatrix(state->obs, 1, nprot_dbase, 1, MAXRES);
    free_cmatrix(state->title_obs, 1, nprot_dbase, 1, MAXLINE);
    free_cmatrix(state->title_seq, 1, nprot_dbase, 1, MAXLINE);
    free_ivector(state->temp, 1, nprot_dbase);
    free_ivector(state->nres, 1, nprot_dbase);

    if (state->resultLength){
        free_cvector(state->predi, 1, state->resultLength);
        free_matrix(state->probai, 1, state->resultLength, 1, 3);
    }

    free(state);
}

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*


void nerror(char error_text[])
/* Error handler */
{
  fprintf(stderr,"run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


int check_boundaries(int indx, int MXINDX, char name[], int maxv)
/* If maxv = 1 check whether indx is smaller than its maximal allowed value,
   if maxv = 0 check whether indx is greater than its minimal allowed value */
{
  if(maxv) {
    if(indx > MXINDX) {
      fprintf(stderr,"Warning: the upper boundary has been reached. Increase %s= %d\n",name,indx);
      return(1);
    }
  } else {
    if(indx < MXINDX) {
      fprintf(stderr,"Warning: the lower boundary has been reached. Decrease %s= %d\n",name,indx);
      return(1);
    }
  }
  return(0);
}


float *vector(long nl, long nh)
/* Allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v = (float *) malloc((size_t) ((nh-nl+1+NR_END) * sizeof(float)));
  if(!v) nerror("allocation failure in vector()");
  return v-nl+NR_END;
}


int *ivector(long nl, long nh)
/* Allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v = (int *) malloc((size_t) ((nh-nl+1+NR_END) * sizeof(int)));
  if(!v) nerror("allocation failure in ivector()");
  return v-nl+NR_END;
}


char *cvector(long nl, long nh)
/* Allocate a char vector with subscript range v[nl..nh] */
{
  char *v;

  v = (char *) malloc((size_t) ((nh-nl+1+NR_END) * sizeof(char)));
  if(!v) nerror("allocation failure in cvector()");
  return v-nl+NR_END;
}


long *lvector(long nl, long nh)
/* Allocate an long vector with subscript range v[nl..nh] */
{
  long *v;

  v = (long *) malloc((size_t) ((nh-nl+1+NR_END) * sizeof(long)));
  if(!v) nerror("allocation failure in lvector()");
  return v-nl+NR_END;
}

unsigned short *svector(long nl, long nh)
/* Allocate an unsigned short vector with subscript range v[nl..nh] */
{
  unsigned short *v;

  v = (unsigned short *) malloc((size_t) ((nh-nl+1+NR_END) * sizeof(unsigned short)));
  if(!v) nerror("allocation failure in svector()");
  return v-nl+NR_END;
}


double *dvector(long nl, long nh)
/* Allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v = (double *) malloc((size_t) ((nh-nl+1+NR_END) * sizeof(double)));
  if(!v) nerror("allocation failure in dvector()");
  return v-nl+NR_END;
}


float **matrix(long nrl, long nrh, long ncl, long nch)
/* Allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{

  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **) malloc((size_t) ((nrow+NR_END) * sizeof(float *)));
  if(!m) nerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *) malloc((size_t) ((nrow*ncol+NR_END) * sizeof(float)));
  if(!m[nrl]) nerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++)
    m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{

  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc((size_t) ((nrow+NR_END) * sizeof(double *)));
  if(!m) nerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *) malloc((size_t) ((nrow*ncol+NR_END) * sizeof(double)));
  if(!m[nrl]) nerror("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++)
    m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}


int **imatrix(long nrl, long nrh, long ncl, long nch)
/* Allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{

  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc((size_t) ((nrow+NR_END) * sizeof(int *)));
  if(!m) nerror("allocation failure 1 in imatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (int *) malloc((size_t) ((nrow*ncol+NR_END) * sizeof(int)));
  if(!m[nrl]) nerror("allocation failure 2 in imatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++)
    m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}


char **cmatrix(long nrl, long nrh, long ncl, long nch)
/* Allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
{

  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  char **m;

  /* allocate pointers to rows */
  m = (char **) malloc((size_t) ((nrow+NR_END) * sizeof(char *)));
  if(!m) nerror("allocation failure 1 in cmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (char *) malloc((size_t) ((nrow*ncol+NR_END) * sizeof(char)));
  if(!m[nrl]) nerror("allocation failure 2 in cmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++)
    m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}


float **submatrix(float **a, long oldrl, long oldrh, long oldcl,  long oldch, long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
  float **m;

  /* allocate array of pointers to rows */
  m = (float **) malloc((size_t) ((nrow+NR_END) * sizeof(float*)));
  if(!m) nerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i = oldrl, j = newrl; i <= oldrh; i++, j++)
    m[j] = a[i] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}


float **convert_matrix(float *a, long nrl, long nrh, long ncl,  long nch)
/* Allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
   declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
   and ncol=nch-ncl+1. The routine should be called with the address &a[0][0]
   as the first argument. */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl;
  float **m;

  /* allocate array of pointers to rows */
  m = (float **) malloc((size_t) ((nrow+NR_END) * sizeof(float*)));
  if(!m) nerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl] = a - ncl;
  for(i = 1, j = nrl + 1; i < nrow; i++, j++)
    m[j] = m[j-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}


float ***f3tensor(long nrl, long nrh, long ncl,  long nch, long ndl, long ndh)
/* Allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  float ***t;

  /* allocate pointers to  pointers to rows */
  t = (float ***) malloc((size_t) ((nrow+NR_END) * sizeof(float**)));
  if(!t) nerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = (float **) malloc((size_t) ((nrow*ncol+NR_END) * sizeof(float *)));
  if(!t[nrl]) nerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = (float *) malloc((size_t) ((nrow*ncol*ndep+NR_END) * sizeof(float)));
  if(!t[nrl][ncl]) nerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for(i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol * ndep;
    for(j = ncl + 1; j <= nch; j++) 
      t[i][j] = t[i][j-1] + ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}


void free_vector(float *v, long nl, long nh)
/* Free a float vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_ivector(int *v, long nl, long nh)
/* Free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_cvector(char *v, long nl, long nh)
/* Free a character vector allocated with cvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_lvector(long *v, long nl, long nh)
/* Free a long vector allocated with lvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_svector(unsigned short *v, long nl, long nh)
/* Free an unsigned short vector allocated with svector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_dvector(double *v, long nl, long nh)
/* Free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* Free a float matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* Free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* Free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* Free a char matrix allocated by cmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* Free a submatrix allocated by submatrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}


void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* Free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}


void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* Free a float 3tensor allocated by f3rensor */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}


int INDMAXVAL(float val[], int i1, int i2)
{
/*
 * Return the index of the greatest element of val between positions i1 and i2
 */
  int i, ini;

  ini = i1;
  for(i = i1+1; i <= i2; i++) {
    if(val[ini] < val[i])
      ini = i;
  }
  return(ini);
}


int INDMINVAL(float val[], int i1, int i2)
{
/*
 * Return the index of the smallest element of val between positions i1 and i2
 */
  int i, ini;

  printf("val[0]= %f\n",val[0]);
  ini = i1;
  for(i = i1+1; i <= i2; i++) {
    if(val[ini] > val[i])
      ini = i;
  }
  return(ini);
}


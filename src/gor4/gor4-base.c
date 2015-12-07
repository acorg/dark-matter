#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>

#ifndef _CFFI_
#include "defines.h"
#include "nrutil.h"
#include "gor4-base.h"
#endif

char buffer[BUFSIZE];

int seq_indx(int c);
int obs_indx(int c);
int read_file(char *fname, int nprot, char **obs, char **title, int *pnter);
void read_fasta(FILE *fp, char **title, char **seq, int *nres, int nprot);
void Parameters(int nprot_dbase, int *nres, char **obs, char **seq);
void predic(int nres, char *seq, char *pred, float **proba);
void First_Pass(int nres, float **proba, char *pred);
void Second_Pass(int nres, float **proba, char *pred);
void printout(int nres, char *seq, char *predi, char *title, float **proba, FILE *fp);



/*
 * External variables
 */

char conf[5] = {' ','H','E','C','S'};
double infopair[3][NPAIRS+1][23][23];
double infodir[3][WINSIZ+1][23];
float ExpInv = 1.0 / (float) WINSIZ;
float nS[4], pS[4];

int main(int argc, char *argv[])
{

/***************************************************************************/
/*                                                                         */
/* GOR secondary structure prediction method version IV                    */
/* J. Garnier, J.-F. Gibrat, B. Robson, Methods in Enzymology,             */
/* R.F. Doolittle Ed., vol 266, 540-553, (1996)                            */
/*                                                                         */
/* For any information please contact: J. Garnier or J.-F. Gibrat          */
/* Unite de Bio-informatique, Batiment des Biotechnologies, I.N.R.A.,      */
/* 78351, Jouy-en-Josas, FRANCE.                                           */
/* tel: 33 (1) 34 65 25 67                                                 */
/* fax: 33 (1) 34 65 22 73                                                 */
/* Email: gibrat@proline.jouy.inra.fr                                      */
/*                                                                         */
/* This program gets its input from the command line                       */
/*                                                                         */
/* This program provides 64.4% residues correctly predicted for a 3 states */
/* prediction (H, E, C) for the 267 proteins of the database using a jack- */
/* knife.                                                                  */
/* Last modification: 07/10/97                                             */
/*                                                                         */
/***************************************************************************/

  FILE *fp, *fp2;
  char **obs, **seq, **SEQ;
  char **title_obs, **title_seq, **TITLE;
  char *Fname1, *Fname2, *Fname3, *Fname4;
  int nprot, nprot_dbase;
  int *temp, *nres, *NRES;
  int i;
  int pro;
  int nerr;
  char *predi;
  float **probai;

  Fname1 = cvector(0,200);
  Fname2 = cvector(0,200);
  Fname3 = cvector(0,200);
  Fname4 = cvector(0,200);

/*
 * Default values for the various parameters
 */

  strcpy(Fname1,"No input value yet");
  strcpy(Fname2,"DATABASE/New_KS.267.seq");
  strcpy(Fname3,"DATABASE/New_KS.267.obs");
  strcpy(Fname4,"No input value yet");

/*
 * Get parameter values from the command line
 */
  
  if(argc == 1) {
    printf("\nUsage %s -prd Fname1 [-seq Fname2] [-obs Fname3] [-pro Fname4]\n",argv[0]);
    printf("Fname1: name of the file containing sequence[s] to be predicted (mandatory)\n");
    printf("Fname2: name of the file containing Kabsch-Sander sequence database (%s)\n",Fname2);
    printf("Fname3: name of the file containing Kabsch-Sander observed secondary structure database (%s)\n",Fname3);
    printf("Fname4: name of the output file that will contain GOR probabilities\n\n\n");
    exit(2);
  }

  if((argc-1)% 2 != 0) {
    printf("Each argument must be preceded by a tag, e.g., -seq foobar.seq\n");
    exit(1);
  }

  for(i = 1; i < argc; i += 2) {
    if(argv[i][0] != '-') {
      printf("invalid tag %s\n",argv[i]);
      exit(1);
    }
  }
  
  for(i = 1; i < argc; i += 2) {
    if(strcmp(argv[i],"-prd") == 0) {
      strcpy(Fname1,argv[i+1]);
    } else if(strcmp(argv[i],"-seq") == 0) {
      strcpy(Fname2,argv[i+1]);
    } else if(strcmp(argv[i],"-obs") == 0) {
      strcpy(Fname3,argv[i+1]);
    } else if(strcmp(argv[i],"-pro") == 0) {
      strcpy(Fname4,argv[i+1]);
    } else {
      printf("Unknown tag: %s\n",argv[i]);
      exit(1);
    }
  }

  if(strcmp(Fname1,"No input value yet") == 0) {
    printf("The name of a file containing sequence[s] to be predicted is mandatory\n");
    exit(1);
  }

  if(strcmp(Fname4,"No input value yet") == 0) {
    fp2 = NULL;
  } else {
    if((fp2 = fopen(Fname4,"w")) == NULL) {
      printf("Unable to open file %s\n",Fname4);
      exit(1);
    }
  }

/*
 * Determine the number of proteins in the Kabsch-Sander data base
 */

  if((fp = fopen(Fname2,"r")) == NULL) {
    printf("Unable to open file %s\n",Fname2);
    exit(1);
  }
  
  nprot_dbase = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '>' || buffer[0] == '!') nprot_dbase++;
  }
  rewind(fp);
  fclose(fp);
  printf("There are %d proteins in Kabsch-Sander database\n\n",nprot_dbase);

/*
 * Memory allocations
 */

  seq = cmatrix(1,nprot_dbase,1,MAXRES);
  obs = cmatrix(1,nprot_dbase,1,MAXRES);
  title_obs = cmatrix(1,nprot_dbase,1,MAXLINE);
  title_seq = cmatrix(1,nprot_dbase,1,MAXLINE);
  temp = ivector(1,nprot_dbase);
  nres = ivector(1,nprot_dbase);
  predi = cvector(1,MAXRES);
  probai = matrix(1,MAXRES,1,3);

/*
 * Input the sequences and observed secondary structures for the data base
 */

  if (!read_file(Fname2,nprot_dbase,seq,title_seq,temp) ||
      !read_file(Fname3,nprot_dbase,obs,title_obs,nres)){
    exit(1);
  }

/*
 * Check that the data are consistent in the two files
 */

  nerr = 0;
  for(i = 1; i <= nprot_dbase; i++)
    if(temp[i] != nres[i]) {
      printf("%dth protein temp= %d nres= %d\n",i,temp[i],nres[i]);
      printf("%s\n%s\n\n",title_seq[i],title_obs[i]);
      nerr++;
    }
  
  for(i = 1; i <= nprot_dbase; i++)
    if(strncmp(title_seq[i],title_obs[i],100) != 0) {
      printf("\n%dth data base protein\n %s \n %s \n",i,title_seq[i],title_obs[i]);
      nerr++;
    }

  if(nerr > 0) {
    printf("%d errors\n",nerr);
    exit(1);
  }

/*
 * Input the sequences for the protein to be predicted
 */

  if((fp = fopen(Fname1,"r")) == NULL) {
    printf("Unable to open file %s\n",Fname1);
    exit(1);
  }

  nprot = 0;
  while(fgets(buffer,BUFSIZE,fp) != NULL) {
    if(buffer[0] == '>' || buffer[0] == '!') nprot++;
  }

  NRES = ivector(1,nprot);
  TITLE = (char **) malloc((size_t) (nprot+1)*sizeof(char *));
  SEQ   = (char **) malloc((size_t) (nprot+1)*sizeof(char *));
  
  rewind(fp);

  read_fasta(fp,TITLE,SEQ,NRES,nprot);

  fclose(fp);

/*
 * Calculate the parameters
 */

  Parameters(nprot_dbase,nres,obs,seq);

/*
 * Predict the secondary structure of protein pro.
 */

  for(pro = 1; pro <= nprot; pro++) {

/*
 * Carry out the prediction for the sequence alone
 */

    predic(NRES[pro],SEQ[pro],predi,probai);
    First_Pass(NRES[pro],probai,predi);
    Second_Pass(NRES[pro],probai,predi); 

/*
 * Print the results for the protein
 */

    printout(NRES[pro],SEQ[pro],predi,TITLE[pro],probai,fp2); 

    if(pro % 100 == 0) printf("%d proteins have been predicted so far\n",pro);
      
  } 

/*
 * Free memory
 */

  free_cmatrix(seq,1,nprot_dbase,1,MAXRES);
  free_cmatrix(obs,1,nprot_dbase,1,MAXRES);
  free_cmatrix(title_obs,1,nprot_dbase,1,MAXLINE);
  free_cmatrix(title_seq,1,nprot_dbase,1,MAXLINE);
  free_ivector(temp,1,nprot_dbase);
  free_ivector(nres,1,nprot_dbase);
  free_cvector(predi,1,MAXRES);
  free_matrix(probai,1,MAXRES,1,3);

  if(fp2 != NULL) {
    fclose(fp2);
  }

  return(0);

}

/*****************************************************************************/
/*                                                                           */
/* This function returns an integer for each amino acid type.                */
/*                                                                           */
/*****************************************************************************/
int seq_indx(int c)
{

  switch(c) {
  case 'A': return(1);
  case 'C': return(2);
  case 'D': return(3);
  case 'E': return(4);
  case 'F': return(5);
  case 'G': return(6);
  case 'H': return(7);
  case 'I': return(8);
  case 'K': return(9);
  case 'L': return(10);
  case 'M': return(11);
  case 'N': return(12);
  case 'P': return(13);
  case 'Q': return(14);
  case 'R': return(15);
  case 'S': return(16);
  case 'T': return(17);
  case 'V': return(18);
  case 'W': return(19);
  case 'Y': return(20);
  case '^': return(21);
  case '-': return(22);
  default : return(23);
  }

}
/*****************************************************************************/
/*                                                                           */
/* This function returns an integer for each secondary structure type.       */
/*                                                                           */
/*****************************************************************************/
int obs_indx(int c)
{

  switch(c) {
  case 'H': return(1);
  case 'E': return(2);
  case 'C': return(3);
  case 'X': return(0);
  }

  return -1;
}
/***********************************************************************************************************/
/*                                                                                                         */
/* This routine performs the prediction of the current protein                                             */
/*                                                                                                         */
/***********************************************************************************************************/
void Normalize(float proba[], double v[]);

void predic(int nres, char *seq, char *pred, float **proba)
{
  double it[3];
  int aa1, aa2;
  int konf, ires;
  int dis1, dis2, np;

/*
 * Calculate sum of information values for each secondary structure type (konf)
 */

  for(ires = 1; ires <= nres; ires++) {  
    it[1] = it[2] = 0.0;
    for(dis1 = -DISLOCATION; dis1 <= +DISLOCATION; dis1++) {
      if(ires+dis1 < 1 || ires +dis1 > nres) {
        if(SKIP_BLANK) continue;                 /* If SKIP_BLANK "amino acid" of type ' ', i.e., */
        aa1 = 21;                                /* aa1 = 21 are not included in the calculation  */
      } else {
        aa1 = seq_indx(seq[ires+dis1]);
      }
      for(dis2 = dis1+1; dis2 <= +DISLOCATION; dis2++) {
	if(ires+dis2 < 1 || ires +dis2 > nres) {
	  if(SKIP_BLANK) continue;                 
	  aa2 = 21;                                
	} else {
	  aa2 = seq_indx(seq[ires+dis2]);
	}
	np = (dis1+8) * (WINSIZ-1) - ((dis1+8)*(dis1+9)/2) + (dis2+8);
	for(konf = 1; konf <= 2; konf++) 
	  it[konf] = it[konf] + infopair[konf][np][aa1][aa2];
      }
    }
    for(dis1 = -DISLOCATION; dis1 <= +DISLOCATION; dis1++) {
      if(ires+dis1 < 1 || ires +dis1 > nres) {
        if(SKIP_BLANK) continue;                 
        aa1 = 21;                                
      } else {
        aa1 = seq_indx(seq[ires+dis1]);
      }
      for(konf = 1; konf <= 2; konf++) {
	it[konf] = it[konf] + infodir[konf][dis1+9][aa1];
      }
    }

    Normalize(proba[ires],it);
    pred[ires] = conf[INDMAXVAL(proba[ires],1,3)];

  }

/*
 * If "blank residues" are not included the first Nterm and the last Cterm residues are predicted as coils
 */
    
  if(SKIP_BLANK) {                              
    for(ires = 1; ires <= Nterm; ires++)        
      pred[ires] = 'C';
    for(ires = nres-Cterm+1; ires <= nres; ires++)
      pred[ires] = 'C';
  }

}
/***************************************************************************/
/*                                                                         */
/*                           Routine Parameters                            */
/*                                                                         */
/***************************************************************************/
int seq_indx(int c);
int obs_indx(int c);
void Indices(int np, int *dis1, int *dis2);

void Parameters(int nprot_dbase, int *nres, char **obs, char **seq)
{
/*
 * Compute the frequencies from proteins in the data base.
 *
 */

  int pro;
  int ires;
  int konf, dis, aa1, aa2, np;
  int dis1, dis2;
  float Singlet[4][WINSIZ+1][23];
  float Doublet[4][NPAIRS+1][23][23];
  double C1, C2;
  float f1, f2, f3;

  C1 = 2. * ExpInv;
  C2 = 1. - C1;

/*
 * Initialisation
 */

  for(konf = 0; konf < 4; konf++)
    for(dis = 0; dis < WINSIZ+1; dis++)
      for(aa1 = 0; aa1 < 23; aa1++) 
	Singlet[konf][dis][aa1] = 0.0;


  for(konf = 0; konf < 4; konf++)
    for(np = 0; np < NPAIRS+1; np++) 
      for(aa1 = 0; aa1 < 23; aa1++)
	for(aa2 = 0; aa2 < 23; aa2++)
	  Doublet[konf][np][aa1][aa2] = 0.0;

  nS[0] = nS[1] = nS[2] = nS[3] = 0;
  
/*
 * Loop over all the proteins of the data base. 
 */

  for(pro = 1; pro <= nprot_dbase; pro++) {

/*
 * Determine frequencies related to the sequence of the query protein (the 1st row in the alignment)
 */

    for(ires = 1; ires <= nres[pro]; ires++) {

      konf = obs_indx(obs[pro][ires]);
      if(konf == 0)                    /* Skip X conformations, i.e., residues for */
	continue;                      /* which the secondary structure is unknown */
	  
      nS[konf]++;   

      for(dis = -DISLOCATION; dis <= DISLOCATION; dis++) {
	if(ires+dis < 1 || ires+dis > nres[pro])
	  aa1 = BLANK;
	else
	  aa1 = seq_indx(seq[pro][ires+dis]);
	Singlet[konf][dis+OFFSET][aa1] += 1.0;
      }
      
      np = 0;
      for(dis1 = -DISLOCATION; dis1 <= DISLOCATION; dis1++) {
	if(ires+dis1 < 1 || ires+dis1 > nres[pro])
	  aa1 = BLANK;
	else
	  aa1 = seq_indx(seq[pro][ires+dis1]);
	for(dis2 = dis1+1; dis2 <= DISLOCATION; dis2++) {
	  if(ires+dis2 < 1 || ires+dis2 > nres[pro])
	    aa2 = BLANK;
	  else
	    aa2 = seq_indx(seq[pro][ires+dis2]);
	  np++;
	  Doublet[konf][np][aa1][aa2] += 1.0;
	}
      }

    }

  } /* End of loop over the proteins in the data base index pro */

/*
 * Calculate probabilities for the 3 secondary structures, H, E and C.
 */

  nS[0] = nS[1] + nS[2] + nS[3];

  for(konf = 1; konf <= 3; konf++)
    pS[konf] = (float) nS[konf] / (float) nS[0];

/*
 * Calculate information parameters (sort of)
 */

  for(konf = 1; konf <= 2; konf++) {
    for(np = 1; np <= NPAIRS; np++) {
      for(aa1 = 1; aa1 <= 21; aa1++) {
	for(aa2 = 1; aa2 <= 21; aa2++) {
	  f1 = Doublet[konf][np][aa1][aa2];
	  f2 = Doublet[3][np][aa1][aa2];
	  if(f1 < MINFREQ) {
	    Indices(np,&dis1,&dis2);
	    f3 = Singlet[konf][dis1][aa1] * Singlet[konf][dis2][aa1] / (float) nS[konf];
	    f1 = (f3 - f1) * interpol_coeff + f1;
	    if(f1 < 1.e-6) f1 = 1.0;
	  }
	  if(f2 < MINFREQ) {
	    Indices(np,&dis1,&dis2);
	    f3 = Singlet[3][dis1][aa1] * Singlet[3][dis2][aa1] / (float) nS[3];
	    f2 = (f3 - f2) * interpol_coeff + f2;
	    if(f2 < 1.e-6) f2 = 1.0;
	  }
	  infopair[konf][np][aa1][aa2] = C1 * (log(f1)-log(f2));
	}
      }
    }
  }

  for(konf = 1; konf <= 2; konf++) {
    for(dis = 1; dis <= WINSIZ; dis++) {
      for(aa1 = 1; aa1 <= 21; aa1++) {
	f1 = Singlet[konf][dis][aa1];
	f2 = Singlet[3][dis][aa1];
	if(f1 < 1.e-6) f1 = 1.0;
	if(f2 < 1.e-6) f2 = 1.0;
	infodir[konf][dis][aa1] = C2 * (log(f2)- log(f1));
      }
    }
  }
  
}

/*****************************************************************************/
/*                                                                           */
/* Determine indices dis1 dis2 as a function of np                           */
/*                                                                           */
/*****************************************************************************/
void Indices(int np, int *dis1, int *dis2)
{
  int i, j, k;

  k = 0;
  for(i = -DISLOCATION; i <= DISLOCATION; i++) {
    for(j= i+1; j <= DISLOCATION; j++) {
      k++;
      if(k == np) {
	*dis1 = i;
	*dis2 = j;
	return;
      }
    }
  }
  printf("Error invalid value of np= %d\n",np);
  exit(1);
}
/*********************************************************************************/
/*                                                                               */
/*                          Normalize the probabilities                          */
/*                                                                               */
/*********************************************************************************/
void Normalize(float proba[], double v[])
{
  double denom;
  
  denom = 1.0 / (1.0 + exp(v[1]) + exp(v[2]));
  proba[1] = exp(v[1]) * denom;
  proba[2] = exp(v[2]) * denom;
  proba[3] = denom;
}
/**********************************************************************************/
/*                                                                                */
/* Print out the results for the current protein.                                 */
/*                                                                                */
/**********************************************************************************/
void printout(int nres, char *seq, char *predi, char *title, float **proba, FILE *fp2)
{
  int ires;
  int nlines, nl;

/*
 * Print the results for the current protein
 */

  printf("\n\n>%s\n",title+1);
  nlines = nres / 50 + 1;

  for(nl = 1; nl < nlines; nl++) {

    for(ires = (nl-1)*50+1; ires <= nl*50; ires++) {
      printf("%c",seq[ires]);
      if(ires % 10 == 0) printf("%c",' ');
    }
    printf("    %s\n","Sequence");

    for(ires = (nl-1)*50+1; ires <= nl*50; ires++) {
      printf("%c",predi[ires]);
      if(ires % 10 == 0) printf("%c",' ');
    }
    printf("    %s\n","Predicted Sec. Struct.");

    printf("\n");

  }

  for(ires = (nlines-1)*50+1; ires <= nlines*50; ires++) { /* last -likely incomplete- line */
    if(ires <= nres) {
      printf("%c",seq[ires]);
    } else {
      printf("%c",' ');
    }
    if(ires % 10 == 0) printf("%c",' ');
  }
  printf("    %s\n","Sequence");

  for(ires = (nlines-1)*50+1; ires <= nlines*50; ires++) {
    if(ires <= nres) {
      printf("%c",predi[ires]);
    } else {
      printf("%c",' ');
    }
    if(ires % 10 == 0) printf("%c",' ');
  }
  printf("    %s\n","Predicted Sec. Struct.");

  printf("\n\n");

  if(fp2 != NULL) {
    fprintf(fp2,"\n\n%s\n%d\n",title+1,nres);
    fprintf(fp2,"SEQ PRD   H     E     C\n");
    for(ires = 1; ires <= nres; ires++)
      fprintf(fp2," %c   %c  %5.3f %5.3f %5.3f\n",seq[ires],predi[ires],proba[ires][1],proba[ires][2],proba[ires][3]);
      
  }

}
/***************************************************************************/
/*                                                                         */
/*                           Routine First_Pass                            */
/*                                                                         */
/***************************************************************************/
int obs_indx(int c);

void First_Pass(int nres, float **proba, char *pred)
{
/*
 * 1) Look for areas that are a mixture of Es and Hs.
 * 2) When such an area is isolated check whether Es and Hs occurs in two blocks.
 * If yes and number of Hs > 4 and number of Es > 3 do nothing
 * In all other cases compute the product of probabilities for all residues in the area
 * and assign to this area the conformation having the highest probability over the area.
 *
 */

  int ires;
  int lim1, lim2;
  int open;
  int kk;
  int type;
  int block[3];
  int nseg;
  int size[3] = {0,4,3};
  double ptot[3];

  pred[1] = pred[nres] = 'C';
  open = 0;
  for(ires = 1; ires <= nres; ires++) {
    if(pred[ires] != 'C') {
      if(!open) {
	open = 1;
	lim1 = ires;
      }
    } else {
      if(open) {
	open = 0;
	lim2 = ires - 1;
	type = obs_indx(pred[lim1]);
	block[1] = block[2] = 0;
	nseg = 1;
	block[nseg]++;
	for(kk = lim1+1; kk <= lim2; kk++) {
	  if(obs_indx(pred[kk]) != type)
	    nseg++;
	  if(nseg <= 2) block[nseg]++;
	  type = obs_indx(pred[kk]);
	}
	if(nseg > 2 || block[1] < size[obs_indx(pred[lim1])] || block[2] < size[obs_indx(pred[lim2])]) {
	  ptot[1] = ptot[2] = 1.0;
	  for(kk = lim1; kk <= lim2; kk++) {
	    ptot[1] = ptot[1] * proba[kk][1];
	    ptot[2] = ptot[2] * proba[kk][2];
	  }
	  if(ptot[1] > ptot[2]) {
	    for(kk = lim1; kk <= lim2; kk++)
	      pred[kk] = 'H';
	  } else {
	    for(kk = lim1; kk <= lim2; kk++)
	      pred[kk] = 'E';
	  }
	}
      }
    }
  }

}
/***************************************************************************/
/*                                                                         */
/*                           Routine Second_Pass                           */
/*                                                                         */
/***************************************************************************/
int obs_indx(int c);

void Second_Pass(int nres, float **proba, char *pred)
{
/*
 * Correct strands having less than 2 and helices having less than 4 residues.
 * Either the secondary structure element is suppressed or additional
 * residues are recruted to reach the required number.
 * 
 */

  int ires, ires1;
  int len;
  int standard[4] = {0,4,2,0};
  int missing;
  int k;
  int lim1, lim2, lim3, lim4, Lim1, Lim2, Lim3, Lim4, KeepNterm, KeepCterm;
  float cost, costmax;
  int type;
  int type_Cterm, type_Nterm;

  len = 0;
  type = obs_indx(pred[1]);
  for(ires = 2; ires <= nres; ires++) {
    if(type != obs_indx(pred[ires])) {
      if(len < standard[type]) { /* Check all possibilities */
	costmax = 0.0;
	missing = standard[type] - len;
/*
 * Check the cost of increasing the secondary structure element
 */
	lim1 = ires - len - missing;
	for(k = 1; k <= missing+1; k++) {
	  lim2 = lim1 + standard[type] - 1;
	  if(lim1 < 1 || lim2 > nres) {
	    lim1++;
	    continue;
	  }
	  cost = 1.0;
	  for(ires1 = lim1; ires1 <= lim2; ires1++)
	    cost *= proba[ires1][type];
	  if(cost > costmax) {
	    costmax = cost;
	    Lim1 = lim1;
	    Lim2 = lim2;
	    KeepNterm = type;
	    Lim3 = 0;
	    Lim4 = -1;
	  }
	  lim1++;
	}
/*
 * Check the cost of suppressing the secondary structure element using the same segments as previously
 */
	type_Nterm = obs_indx(pred[ires-len-1]);
	type_Cterm = obs_indx(pred[ires]);
	lim1 = ires - len - missing;
	for(k = 1; k <= missing+1; k++) {
	  lim4 = lim1 + standard[type] - 1;
	  if(lim1 < 1 || lim4 > nres) {
	    lim1++;
	    continue;
	  }
	  lim2 = ires - 1;
	  lim3 = lim2 + 1;
	  while(lim3 >= ires - len) {
	    cost = 1.0;
	    for(ires1 = lim1; ires1 <= lim2; ires1++)
	      cost *= proba[ires1][type_Nterm];
	    for(ires1 = lim3; ires1 <= lim4; ires1++)
	      cost *= proba[ires][type_Cterm];
	    if(cost > costmax) {
	      costmax = cost;
	      Lim1 = lim1;
	      Lim2 = lim2;
	      Lim3 = lim3;
	      Lim4 = lim4;
	      KeepNterm = type_Nterm;
	      KeepCterm = type_Cterm;
	    }
	    lim2--;
	    lim3--;
	  }
	  lim1++;
	}
/*
 * Modify pred accordingly
 */
	for(ires1 = Lim1; ires1 <= Lim2; ires1++)
	  pred[ires1] = conf[KeepNterm];
	for(ires1 = Lim3; ires1 <= Lim4; ires1++)
	  pred[ires1] = conf[KeepCterm];
/*
 * Move to the end of the modified segment if necessary
 */
	if(Lim2 > ires || Lim4 > ires) {
	  if(Lim2 > Lim4)
	    ires = Lim2;
	  else 
	    ires = Lim4;
	}

      } /* End of segment correction */
 
      len = 1;
    } else {
      len++;
    }
    type = obs_indx(pred[ires]);
  }

}
#include <stdio.h>
#include <string.h>

/***************************************************************************/
/*                                                                         */
/*                           Routine read_fasta                            */
/*                                                                         */
/***************************************************************************/
void read_fasta(FILE *fp, char **title, char **seq, int *nres, int nprot)
{
/*
 * Reads sequences in FASTA format (>title\n sequence) or GOR/HOMOL format (!title\n sequence ending with @)
 * The number of sequences needs not be specified in advance.
 * The one letter-code in uppercase is used for sequences. X is allowed for non conventional amino acid (though it
 * generates a warning). All other characters are ignored.
 */

  int open_title;
  int Tindx;
  int nr;
  int c;
  int i;
  int np;

  open_title = nr = Tindx = np = 0;
  while((c=fgetc(fp)) != EOF) {
    if(c == '!' || c == '>') {               /* Beginning of a new protein: start of the title */
      if(np > 0) {
	nres[np] = nr;
	seq[np] = (char *) malloc((size_t) (nr+2)*sizeof(char));
	for(i = 1; i <= nr; i++)
	  seq[np][i] = buffer[i];
	seq[np][nr+1] = '\0';
      }
      open_title = 1;
      (np)++;
      if(np > nprot) {
	printf("Warning the program is supposed to read %d proteins and so far %d proteins have been read\n",nprot,np);
	exit(1);
      }
    } else if(c == '\n' && open_title) {     /* End of title for a protein */
      open_title = 0;
      title[np] = (char *) malloc((size_t) (Tindx+2)*sizeof(char));
      for(i = 1; i <= Tindx; i++)
	title[np][i] = buffer[i];
      title[np][++Tindx] = '\0';
      Tindx = 0;
      nr = 0;
    } else {
      if(open_title) {                       /* Reads title characters */
	Tindx++;
	if(Tindx > BUFSIZE - 1) {
	  printf("Increase the value of BUFSIZE\n");
	  exit(1);
	}
	buffer[Tindx] = c;
      } else {                               /* Reads protein characters */
	if((c >= 'A' && c < 'Z') && c != 'B' && c != 'J' && c != 'O' && c != 'U') {
	  nr++;
	  if(nr > BUFSIZE) {
	    printf("increase the size of BUFSIZE\n");
	    exit(1);
	  }
	  buffer[nr] = c;
	} else {
	  if(c != '\t' && c != '\n' && c != ' ' && c != '@')
	    printf("Warning non standard amino acid: %c for protein # %d at position %d\n",c,np,nr);
	}
      }
    }
  }
  
  nres[np] = nr;
  seq[np] = (char *) malloc((size_t) (nr+2)*sizeof(char));
  for(i = 1; i <= nr; i++)
    seq[np][i] = buffer[i];
  seq[np][nr+1] = '\0';
  
  printf("%d proteins have been read on sequence file\n\n",np);

}

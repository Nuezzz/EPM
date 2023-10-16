
#ifndef _HAMITONIAN
#define _HAMITONIAN
#include "eigenlapack.h"
#include "atom.h"
#include <complex.h>

double complex *HSO ( Lattice *s, Eigen *d, double *k_vec);
double complex *HLocal  (Lattice *s, Eigen *d);
double complex PotentialMix(Lattice *s, double q[3]);
double complex *H_tot_so ( Lattice *s, Eigen *d, double *Kvec);
double complex  *H_tot_loc   ( Lattice *s, Eigen *d, double *Kvec);
int CalcBand            (Eigen *d, double complex *H, int bands, int N_RANK);


void  PrintEigen        (FILE *fp, double *E, double *k, int N);
FILE  *OpenBandFile     ( char* simname);



void PPtest             (Lattice *s,char* simname);
#endif

#ifndef _HAMITONIAN
#define _HAMITONIAN
#include "eigenlapack.h"
#include "atom.h"
#include <complex.h>


double complex *HLocal  (Lattice *s, Eigen *d);
double complex  *HTot_loc   ( Lattice *s, Eigen *d, double *Kvec);
int CalcBand            (Eigen *d, double complex *H, int bands);

void  PrintEigen        (FILE *fp, double *E, double *k, int N);
FILE  *OpenBandFile     ( char* simname);



void PPtest             (Lattice *s,char* simname);
#endif
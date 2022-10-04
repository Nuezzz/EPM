
#ifndef _HAMITONIAN
#define _HAMITONIAN
#include "eigenlapack.h"
#include "atom.h"
#include <complex.h>


double complex *HLocal(Lattice *s, int NG );
double complex  *HTot ( Lattice *s, int NG, double *Kvec);
int CalcBand( Lattice *s, double *Kvec, int NG, int bands);

void  PrintEigen(FILE *fp, double *E, double *k, int N);
FILE  *OpenBandFile( char* simname);



void PPtest(Lattice *s,char* simname);
#endif

#ifndef _HAMITONIAN
#define _HAMITONIAN
#include "eigenlapack.h"
#include "atom.h"
#include <complex.h>

#define BAND_FILENAME     "band_structure"
#define WAVE_FILENAME     "wave_function"

double complex *HSO ( Lattice *s, Eigen *d, double *k_vec);
double complex *HLocal  (Lattice *s, Eigen *d);
double complex PotentialMix(Lattice *s, double q[3]);
double complex *H_tot_so ( Lattice *s, Eigen *d, double *Kvec);
double complex  *H_tot_loc   ( Lattice *s, Eigen *d, double *Kvec);
int CalcBand            (Eigen *d, double complex *H, int bands, int N_RANK);

void  PrintNG           (FILE *fp, int *NG, double *k, int N);
void  PrintEigen        (FILE *fp, double *E, double *k, int N);
FILE  *OpenBandFile     ( char* pathname, unsigned int label);
FILE  *OpenWaveFile       ( char* pathname, unsigned int label);



void PPtest             (Lattice *s,char* simname);
#endif